//! Rust inner loop for npl's ACT guided-exchange ordering search.
//!
//! Mirrors, step for step, the Python in
//! `npl/optimization/local_optimization/act_exchange_operator.py`,
//! `npl/optimization/local_optimization/local_optimization.py` and
//! `npl/optimization/basin_hopping/basin_hopping.py` for the case that
//! dominates production runtime: `model='ACT'` with a linear (BRR/TEC) energy
//! model on a two-element particle.
//!
//! Everything that loop touches is integer bookkeeping plus one dot product,
//! so no Python object survives the port:
//!
//!   feature(i) = 91 * symbol(i) + cn_offset[cn(i)] + (# neighbours of type a)
//!   E          = sum_i coef[feature(i)]        (the TEC feature vector is a
//!                                               histogram, so the dot product
//!                                               is a sum over atoms)
//!
//! Symbols are 0 for the alphabetically first element and 1 for the second,
//! matching `TopologicalEnvironmentClassifier`'s `symbols.index(symbol)`.

use pyo3::prelude::*;

/// Sorted-by-key list of (exchange energy, atom index).
///
/// `SortedKeyList` in the Python is used for two things only: take the
/// smallest, and index by rank. A sorted `Vec` does both, and its O(n) insert
/// is a memmove of a few kB - far cheaper than the B-tree walk it replaces.
// ponytail: sorted Vec, O(n) insert; swap for an order-statistic tree only if
// particles ever grow past ~10^4 atoms.
struct SortedList {
    items: Vec<(f64, u32)>,
}

impl SortedList {
    fn new() -> Self {
        SortedList { items: Vec::new() }
    }

    fn cmp_key(a: &(f64, u32), b: &(f64, u32)) -> std::cmp::Ordering {
        a.0.total_cmp(&b.0).then(a.1.cmp(&b.1))
    }

    fn insert(&mut self, key: f64, index: u32) {
        let entry = (key, index);
        let at = match self.items.binary_search_by(|p| Self::cmp_key(p, &entry)) {
            Ok(i) | Err(i) => i,
        };
        self.items.insert(at, entry);
    }

    fn remove(&mut self, key: f64, index: u32) {
        let entry = (key, index);
        if let Ok(i) = self.items.binary_search_by(|p| Self::cmp_key(p, &entry)) {
            self.items.remove(i);
        }
    }

    fn len(&self) -> usize {
        self.items.len()
    }

    fn at(&self, rank: usize) -> u32 {
        self.items[rank].1
    }
}

/// xoshiro256++ - a few lines, no crate, and good enough for picking swap
/// candidates. The Python used `random.sample`, so streams differ between the
/// two implementations by construction; nothing downstream depends on
/// reproducing Python's stream, only on the search being unbiased.
struct Rng {
    s: [u64; 4],
}

impl Rng {
    fn new(seed: u64) -> Self {
        // SplitMix64 to spread a single seed over the state
        let mut z = seed;
        let mut s = [0u64; 4];
        for slot in s.iter_mut() {
            z = z.wrapping_add(0x9E3779B97F4A7C15);
            let mut x = z;
            x = (x ^ (x >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
            x = (x ^ (x >> 27)).wrapping_mul(0x94D049BB133111EB);
            *slot = x ^ (x >> 31);
        }
        Rng { s }
    }

    fn next_u64(&mut self) -> u64 {
        let result = self.s[0]
            .wrapping_add(self.s[3])
            .rotate_left(23)
            .wrapping_add(self.s[0]);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];
        self.s[2] ^= t;
        self.s[3] = self.s[3].rotate_left(45);
        result
    }

    /// Uniform in [0, bound), debiased by rejection.
    fn below(&mut self, bound: u64) -> u64 {
        let zone = u64::MAX - u64::MAX % bound;
        loop {
            let v = self.next_u64();
            if v < zone {
                return v % bound;
            }
        }
    }

    fn permutation(&mut self, n: usize) -> Vec<u32> {
        let mut out: Vec<u32> = (0..n as u32).collect();
        for i in (1..n).rev() {
            let j = self.below(i as u64 + 1) as usize;
            out.swap(i, j);
        }
        out
    }
}

/// One particle plus the ACT exchange bookkeeping, all of it flat arrays.
struct Act {
    n: usize,
    offsets: Vec<u32>,
    neighbours: Vec<u32>,
    symbol: Vec<u8>,
    feature: Vec<u32>,
    cn_offset: [u32; 13],
    coef: Vec<f64>,
    /// environment_energies[i] - environment_energies[i + 91], i in 0..91
    env_diff: Vec<f64>,
    /// environment_energies[i+1] - environment_energies[i], with -1000 appended
    neigh_diff: Vec<f64>,
    /// exchange energy currently stored for each atom (the sort key in use)
    key: Vec<f64>,
    s1: SortedList,
    s2: SortedList,
    energy: f64,
}

impl Act {
    fn new(
        offsets: Vec<u32>,
        neighbours: Vec<u32>,
        symbol: Vec<u8>,
        coef: Vec<f64>,
        environment_energies: Vec<f64>,
    ) -> Self {
        let n = symbol.len();
        let n_envs = environment_energies.len() / 2;
        let env_diff: Vec<f64> = (0..n_envs)
            .map(|i| environment_energies[i] - environment_energies[i + n_envs])
            .collect();
        let mut neigh_diff: Vec<f64> = (0..(n_envs * 2 - 1))
            .map(|i| environment_energies[i + 1] - environment_energies[i])
            .collect();
        neigh_diff.push(-1000.0);

        let mut cn_offset = [0u32; 13];
        for cn in 0..13u32 {
            cn_offset[cn as usize] = cn * (cn + 1) / 2;
        }

        let mut act = Act {
            n,
            offsets,
            neighbours,
            symbol,
            feature: vec![0; n],
            cn_offset,
            coef,
            env_diff,
            neigh_diff,
            key: vec![0.0; n],
            s1: SortedList::new(),
            s2: SortedList::new(),
            energy: 0.0,
        };
        act.recompute_all_features();
        act.bind();
        act
    }

    #[inline]
    fn neighbours_of(&self, i: usize) -> &[u32] {
        let lo = self.offsets[i] as usize;
        let hi = self.offsets[i + 1] as usize;
        &self.neighbours[lo..hi]
    }

    #[inline]
    fn compute_feature(&self, i: usize) -> u32 {
        let lo = self.offsets[i] as usize;
        let hi = self.offsets[i + 1] as usize;
        let cn = hi - lo;
        let mut n_a = 0u32;
        for k in lo..hi {
            if self.symbol[self.neighbours[k] as usize] == 0 {
                n_a += 1;
            }
        }
        91 * self.symbol[i] as u32 + self.cn_offset[cn] + n_a
    }

    fn recompute_all_features(&mut self) {
        for i in 0..self.n {
            self.feature[i] = self.compute_feature(i);
        }
    }

    /// Total energy from scratch: the TEC feature vector is a histogram, so the
    /// dot product with the ridge coefficients is a sum over atoms. Recomputed
    /// every step, exactly as the Python does, so no incremental drift.
    fn total_energy(&self) -> f64 {
        let mut e = 0.0;
        for i in 0..self.n {
            e += self.coef[self.feature[i] as usize];
        }
        e
    }

    /// Exchange energy of one atom, given current features everywhere.
    ///
    /// The `neigh_feature - 1` lookup is a signed index in the Python and wraps
    /// to the final -1000 sentinel when the neighbour feature is 0; that wrap is
    /// reproduced here rather than clamped.
    #[inline]
    fn exchange_energy(&self, i: usize) -> f64 {
        let last = self.neigh_diff.len() - 1;
        let mut change = 0.0;
        if self.symbol[i] == 0 {
            for &nb in self.neighbours_of(i) {
                let nf = self.feature[nb as usize] as usize;
                let idx = if nf == 0 { last } else { nf - 1 };
                change -= self.neigh_diff[idx];
            }
            -self.env_diff[(self.feature[i] % self.env_diff.len() as u32) as usize] + change
        } else {
            for &nb in self.neighbours_of(i) {
                change += self.neigh_diff[self.feature[nb as usize] as usize];
            }
            self.env_diff[(self.feature[i] % self.env_diff.len() as u32) as usize] + change
        }
    }

    fn bind(&mut self) {
        self.s1 = SortedList::new();
        self.s2 = SortedList::new();
        for i in 0..self.n {
            let e = self.exchange_energy(i);
            self.key[i] = e;
            if self.symbol[i] == 0 {
                self.s1.insert(e, i as u32);
            } else {
                self.s2.insert(e, i as u32);
            }
        }
        self.energy = self.total_energy();
    }

    /// The two atoms plus every atom whose feature or exchange energy can move:
    /// one shell for the features, a second shell because an exchange energy
    /// reads its neighbours' features.
    fn neighbourhood(&self, i1: usize, i2: usize) -> Vec<u32> {
        let mut seen = vec![false; self.n];
        let mut inner: Vec<u32> = Vec::with_capacity(32);
        let push = |idx: u32, seen: &mut Vec<bool>, out: &mut Vec<u32>| {
            if !seen[idx as usize] {
                seen[idx as usize] = true;
                out.push(idx);
            }
        };
        push(i1 as u32, &mut seen, &mut inner);
        push(i2 as u32, &mut seen, &mut inner);
        for &nb in self.neighbours_of(i1) {
            push(nb, &mut seen, &mut inner);
        }
        for &nb in self.neighbours_of(i2) {
            push(nb, &mut seen, &mut inner);
        }

        let mut all = inner.clone();
        for &x in &inner {
            for &nb in self.neighbours_of(x as usize) {
                push(nb, &mut seen, &mut all);
            }
        }
        all
    }

    /// Swap two atoms and bring features, exchange energies, sort order and the
    /// total energy back in sync.
    fn apply_swap(&mut self, i1: usize, i2: usize) {
        self.symbol.swap(i1, i2);
        let nbhd = self.neighbourhood(i1, i2);

        // drop the stale sort entries first, while `key` still holds the value
        // they were filed under
        for &idx in &nbhd {
            let i = idx as usize;
            if self.symbol[i] == 0 {
                // exchanged atoms already carry their new symbol, so they come
                // out of the list belonging to the old one
                if i == i1 || i == i2 {
                    self.s2.remove(self.key[i], idx);
                } else {
                    self.s1.remove(self.key[i], idx);
                }
            } else if i == i1 || i == i2 {
                self.s1.remove(self.key[i], idx);
            } else {
                self.s2.remove(self.key[i], idx);
            }
        }

        for &idx in &nbhd {
            self.feature[idx as usize] = self.compute_feature(idx as usize);
        }

        for &idx in &nbhd {
            let i = idx as usize;
            let e = self.exchange_energy(i);
            self.key[i] = e;
            if self.symbol[i] == 0 {
                self.s1.insert(e, idx);
            } else {
                self.s2.insert(e, idx);
            }
        }

        self.energy = self.total_energy();
    }

    /// Returns the two swapped atoms and the exchange energy the operator
    /// predicted for the swap, read before the swap as the Python does.
    fn guided_exchange(&mut self) -> (usize, usize, f64) {
        let i1 = self.s1.at(0) as usize;
        let i2 = self.s2.at(0) as usize;
        let flip = self.key[i1] + self.key[i2];
        self.apply_swap(i1, i2);
        (i1, i2, flip)
    }

    fn basin_hop_step(&mut self, rng: &mut Rng) -> (usize, usize) {
        let n1 = self.s1.len();
        let n2 = self.s2.len();
        let wide = n1.max(n2);
        let narrow = n1.min(n2);
        let order = rng.permutation(wide);

        let mut index: isize = -1;
        while index < narrow as isize {
            index += 1;
            // ponytail: the Python indexes `random_indices[index]` here with
            // index == narrow on the fall-through, which is out of range when
            // n1 == n2; clamping keeps equal compositions working and is
            // identical wherever the Python does not raise.
            let pick = order[(index as usize).min(wide - 1)] as usize;
            let i1 = self.s1.at(pick % n1) as usize;
            let i2 = self.s2.at(pick % n2) as usize;
            if self.key[i1] + self.key[i2] > 0.0 {
                self.apply_swap(i1, i2);
                return (i1, i2);
            }
        }

        let pick = order[(index as usize).min(wide - 1)] as usize;
        let i1 = self.s1.at(pick % n1) as usize;
        let i2 = self.s2.at(pick % n2) as usize;
        self.apply_swap(i1, i2);
        (i1, i2)
    }

    /// Largest disagreement between the incrementally maintained exchange
    /// energies and a from-scratch recompute. Used by the test suite: it checks
    /// the hard part of the port without depending on tie-breaking order.
    fn bookkeeping_error(&self) -> f64 {
        let mut worst: f64 = 0.0;
        for i in 0..self.n {
            worst = worst.max((self.key[i] - self.exchange_energy(i)).abs());
        }
        worst
    }
}

fn build(
    offsets: Vec<u32>,
    neighbours: Vec<u32>,
    symbol: Vec<u8>,
    coef: Vec<f64>,
    environment_energies: Vec<f64>,
) -> PyResult<Act> {
    if offsets.len() != symbol.len() + 1 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "offsets must have one more entry than there are atoms",
        ));
    }
    if environment_energies.len() % 2 != 0 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "environment_energies must have an even length",
        ));
    }
    if symbol.iter().any(|&s| s > 1) {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "symbols must be 0 or 1 (two elements only)",
        ));
    }
    Ok(Act::new(offsets, neighbours, symbol, coef, environment_energies))
}

/// Energy of one ordering under the linear TEC model.
#[pyfunction]
fn act_energy(
    offsets: Vec<u32>,
    neighbours: Vec<u32>,
    symbols: Vec<u8>,
    coef: Vec<f64>,
) -> PyResult<f64> {
    let act = build(offsets, neighbours, symbols, coef, vec![0.0; 182])?;
    Ok(act.total_energy())
}

/// Deterministic guided-exchange descent to a local minimum.
///
/// Returns (symbols, accepted energies as (energy, step), bookkeeping error).
#[pyfunction]
fn act_local_optimization(
    offsets: Vec<u32>,
    neighbours: Vec<u32>,
    symbols: Vec<u8>,
    coef: Vec<f64>,
    environment_energies: Vec<f64>,
) -> PyResult<(Vec<u8>, Vec<(f64, u64)>, f64)> {
    let mut act = build(offsets, neighbours, symbols, coef, environment_energies)?;

    let mut start_energy = act.energy;
    let mut accepted: Vec<(f64, u64)> = vec![(start_energy, 0)];
    let mut step: u64 = 0;

    loop {
        step += 2;
        let (i1, i2, _flip) = act.guided_exchange();
        let new_energy = act.energy;
        if new_energy < start_energy {
            start_energy = new_energy;
            accepted.push((new_energy, step));
        } else {
            accepted.push((start_energy, step));
            act.apply_swap(i1, i2);
            break;
        }
    }

    let error = act.bookkeeping_error();
    Ok((act.symbol.clone(), accepted, error))
}

/// ACT guided-exchange basin hopping.
///
/// Returns (best symbols, best energy, lowest energies as (energy, step),
/// flip list as (predicted flip, realised dE, [i1, i2]), bookkeeping error).
#[pyfunction]
#[pyo3(signature = (offsets, neighbours, symbols, coef, environment_energies,
                    attempts, hops, seed))]
fn act_basin_hopping(
    offsets: Vec<u32>,
    neighbours: Vec<u32>,
    symbols: Vec<u8>,
    coef: Vec<f64>,
    environment_energies: Vec<f64>,
    attempts: u64,
    hops: u64,
    seed: u64,
) -> PyResult<(Vec<u8>, f64, Vec<(f64, u64)>, Vec<(f64, f64, [u32; 2])>, f64)> {
    let mut act = build(offsets, neighbours, symbols, coef, environment_energies)?;
    let mut rng = Rng::new(seed);
    let mut flips: Vec<(f64, f64, [u32; 2])> = Vec::new();

    let mut start_energy = act.energy;
    let mut lowest_energy = start_energy;
    let mut best_energy = start_energy;
    let mut best_symbols = act.symbol.clone();
    let mut lowest: Vec<(f64, u64)> = vec![(start_energy, 0)];
    let mut step: u64 = 0;

    for _ in 0..attempts {
        loop {
            step += 1;
            let (i1, i2, flip) = act.guided_exchange();
            let new_energy = act.energy;
            flips.push((flip, new_energy - start_energy, [i1 as u32, i2 as u32]));

            if new_energy < start_energy {
                start_energy = new_energy;
                if new_energy < lowest_energy {
                    lowest_energy = new_energy;
                    lowest.push((lowest_energy, step));
                }
            } else {
                if lowest_energy == start_energy {
                    // the Python snapshots the ordering with this rejected
                    // exchange undone, then puts it back and carries on
                    let mut snapshot = act.symbol.clone();
                    snapshot.swap(i1, i2);
                    best_symbols = snapshot;
                    best_energy = start_energy;
                }
                break;
            }
        }

        for _ in 0..hops {
            step += 1;
            act.basin_hop_step(&mut rng);
            start_energy = act.energy;
        }
    }

    lowest.push((lowest_energy, step));
    let error = act.bookkeeping_error();
    Ok((best_symbols, best_energy, lowest, flips, error))
}


// ----------------------------------------------------------------- genetic
//
// The Python GA spends most of its time not searching: profiling one restart
// on a 711-atom particle showed 20.1 s of `copy.deepcopy` inside the mutation
// operator and 11.4 s rebuilding a neighbour list with ASE inside cut-and-
// splice. Both are pure waste - every member of the population sits on the
// same lattice, so the geometry and the neighbour list never change. Here a
// member is just its symbol vector, and neither cost exists.

/// Energy of an ordering without building the exchange bookkeeping, for
/// scoring population members.
fn energy_of(offsets: &[u32], neighbours: &[u32], symbol: &[u8], coef: &[f64],
             cn_offset: &[u32; 13]) -> f64 {
    let mut e = 0.0;
    for i in 0..symbol.len() {
        let lo = offsets[i] as usize;
        let hi = offsets[i + 1] as usize;
        let mut n_a = 0u32;
        for k in lo..hi {
            if symbol[neighbours[k] as usize] == 0 {
                n_a += 1;
            }
        }
        e += coef[(91 * symbol[i] as u32 + cn_offset[hi - lo] + n_a) as usize];
    }
    e
}

/// Rank-based selection weights for a population sorted best-first: the best
/// member always gets 1 and the worst exp(-pressure), whatever the energies
/// are. Mirrors `rank_fitness`.
fn rank_weights(n: usize, pressure: f64) -> Vec<f64> {
    if n == 1 {
        return vec![1.0];
    }
    (0..n)
        .map(|i| (-pressure * (i as f64 / (n - 1) as f64)).exp())
        .collect()
}

/// One weighted draw, by cumulative sum.
fn weighted_pick(weights: &[f64], rng: &mut Rng) -> usize {
    let total: f64 = weights.iter().sum();
    let mut target = (rng.next_u64() as f64 / u64::MAX as f64) * total;
    for (i, w) in weights.iter().enumerate() {
        target -= w;
        if target <= 0.0 {
            return i;
        }
    }
    weights.len() - 1
}

/// `n` distinct atoms carrying `want`, drawn without replacement.
fn sample_by_symbol(symbol: &[u8], want: u8, n: usize, rng: &mut Rng) -> Vec<usize> {
    let mut pool: Vec<usize> = (0..symbol.len()).filter(|&i| symbol[i] == want).collect();
    let take = n.min(pool.len());
    for i in 0..take {
        let j = i + rng.below((pool.len() - i) as u64) as usize;
        pool.swap(i, j);
    }
    pool.truncate(take);
    pool
}

/// Flip randomly chosen atoms until the ordering has `target_a` atoms of the
/// first element, the two-element case of `adjust_stoichiometry`.
fn fix_stoichiometry(symbol: &mut [u8], target_a: usize, rng: &mut Rng) {
    let count_a = symbol.iter().filter(|&&s| s == 0).count();
    if count_a > target_a {
        for i in sample_by_symbol(symbol, 0, count_a - target_a, rng) {
            symbol[i] = 1;
        }
    } else if count_a < target_a {
        for i in sample_by_symbol(symbol, 1, target_a - count_a, rng) {
            symbol[i] = 0;
        }
    }
}

/// Take the atoms on one side of a random plane from `p1` and the rest from
/// `p2`, then restore `p1`'s composition. Mirrors `CutAndSpliceOperator` with
/// its default random normal direction, minus the particle rebuild: both
/// parents index the same lattice, so the offspring keeps that indexing and
/// inherits the neighbour list unchanged.
fn cut_and_splice(p1: &[u8], p2: &[u8], positions: &[f64], center: [f64; 3],
                  max_radius: f64, rng: &mut Rng) -> Vec<u8> {
    let n = p1.len();
    let unit = |rng: &mut Rng| {
        loop {
            let v = [
                1.0 - 2.0 * (rng.next_u64() as f64 / u64::MAX as f64),
                1.0 - 2.0 * (rng.next_u64() as f64 / u64::MAX as f64),
                1.0 - 2.0 * (rng.next_u64() as f64 / u64::MAX as f64),
            ];
            let norm = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            if norm > 1e-12 {
                return [v[0] / norm, v[1] / norm, v[2] / norm];
            }
        }
    };

    loop {
        let normal = unit(rng);
        let anchor_dir = unit(rng);
        let radius = (rng.next_u64() as f64 / u64::MAX as f64) * max_radius;
        let anchor = [
            anchor_dir[0] * radius + center[0],
            anchor_dir[1] * radius + center[1],
            anchor_dir[2] * radius + center[2],
        ];

        let mut child = vec![0u8; n];
        let (mut n_positive, mut n_negative) = (0usize, 0usize);
        for i in 0..n {
            let d = (positions[3 * i] - anchor[0]) * normal[0]
                + (positions[3 * i + 1] - anchor[1]) * normal[1]
                + (positions[3 * i + 2] - anchor[2]) * normal[2];
            if d > 0.0 {
                child[i] = p1[i];
                n_positive += 1;
            } else {
                // atoms exactly on the plane fall to p2; the Python drops them
                // from both parents, which is the same measure-zero case
                child[i] = p2[i];
                if d < 0.0 {
                    n_negative += 1;
                }
            }
        }
        if n_positive > 0 && n_negative > 0 {
            return child;
        }
    }
}

/// Steady-state genetic algorithm over the orderings of one particle.
///
/// Returns (final population symbols sorted best-first, their energies,
/// best-energy history as (energy, evaluations), evaluations, generations,
/// accepted, duplicates).
#[pyfunction]
#[pyo3(signature = (offsets, neighbours, positions, masses, population, coef,
                    environment_energies, patience, max_offspring_attempts,
                    cut_and_splice_p, selection_pressure, mutation_low,
                    mutation_high, improvement_tol, seed))]
#[allow(clippy::too_many_arguments)]
fn act_genetic_algorithm(
    offsets: Vec<u32>,
    neighbours: Vec<u32>,
    positions: Vec<f64>,
    masses: [f64; 2],
    population: Vec<Vec<u8>>,
    coef: Vec<f64>,
    environment_energies: Vec<f64>,
    patience: u64,
    max_offspring_attempts: u64,
    cut_and_splice_p: f64,
    selection_pressure: f64,
    mutation_low: f64,
    mutation_high: f64,
    improvement_tol: f64,
    seed: u64,
) -> PyResult<(Vec<Vec<u8>>, Vec<f64>, Vec<(f64, u64)>, u64, u64, u64, u64)> {
    if population.is_empty() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "the starting population is empty",
        ));
    }
    let n = population[0].len();
    if positions.len() != 3 * n {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "positions must hold three coordinates per atom",
        ));
    }

    let mut cn_offset = [0u32; 13];
    for cn in 0..13u32 {
        cn_offset[cn as usize] = cn * (cn + 1) / 2;
    }
    let mut rng = Rng::new(seed);

    // geometry is shared by every member, so the centre of mass and the
    // radius are measured once from the first member's composition
    let reference = &population[0];
    let target_a = reference.iter().filter(|&&s| s == 0).count();
    let total_mass: f64 = reference.iter().map(|&s| masses[s as usize]).sum();
    let mut center = [0.0f64; 3];
    for i in 0..n {
        let m = masses[reference[i] as usize];
        for axis in 0..3 {
            center[axis] += m * positions[3 * i + axis];
        }
    }
    for axis in 0..3 {
        center[axis] /= total_mass;
    }
    let mut max_radius: f64 = 0.0;
    for i in 0..n {
        let dx = positions[3 * i] - center[0];
        let dy = positions[3 * i + 1] - center[1];
        let dz = positions[3 * i + 2] - center[2];
        max_radius = max_radius.max((dx * dx + dy * dy + dz * dz).sqrt());
    }

    let descend = |symbols: Vec<u8>| -> (Vec<u8>, f64, u64) {
        let mut act = Act::new(offsets.clone(), neighbours.clone(), symbols,
                               coef.clone(), environment_energies.clone());
        let mut start_energy = act.energy;
        let mut step: u64 = 0;
        loop {
            step += 2;
            let (i1, i2, _flip) = act.guided_exchange();
            if act.energy < start_energy {
                start_energy = act.energy;
            } else {
                act.apply_swap(i1, i2);
                break;
            }
        }
        (act.symbol.clone(), start_energy, step)
    };

    let mut members: Vec<(f64, Vec<u8>)> = population
        .into_iter()
        .map(|s| {
            let e = energy_of(&offsets, &neighbours, &s, &coef, &cn_offset);
            (e, s)
        })
        .collect();
    let mut evaluations = members.len() as u64;
    members.sort_by(|a, b| a.0.total_cmp(&b.0));

    let mut seen: std::collections::HashSet<Vec<u8>> =
        members.iter().map(|m| m.1.clone()).collect();
    let mut history: Vec<(f64, u64)> = vec![(members[0].0, evaluations)];

    let low = ((mutation_low * n as f64).round() as usize).max(1);
    let high = ((mutation_high * n as f64).round() as usize).max(low);

    let (mut generation, mut accepted, mut duplicates) = (0u64, 0u64, 0u64);
    let mut unsuccessful = 0u64;

    while unsuccessful < patience {
        generation += 1;
        let weights = rank_weights(members.len(), selection_pressure);

        let mut offspring: Option<(f64, Vec<u8>)> = None;
        for _ in 0..max_offspring_attempts {
            let candidate = if (rng.next_u64() as f64 / u64::MAX as f64) < cut_and_splice_p
                && members.len() > 1
            {
                let first = weighted_pick(&weights, &mut rng);
                let mut remaining = weights.clone();
                remaining[first] = 0.0;
                let second = weighted_pick(&remaining, &mut rng);
                let mut child = cut_and_splice(&members[first].1, &members[second].1,
                                               &positions, center, max_radius, &mut rng);
                fix_stoichiometry(&mut child, target_a, &mut rng);
                child
            } else {
                let parent = weighted_pick(&weights, &mut rng);
                let mut child = members[parent].1.clone();
                let n_swaps = low + rng.below((high - low + 1) as u64) as usize;
                let firsts = sample_by_symbol(&child, 0, n_swaps, &mut rng);
                let seconds = sample_by_symbol(&child, 1, n_swaps, &mut rng);
                for (a, b) in firsts.into_iter().zip(seconds) {
                    child.swap(a, b);
                }
                child
            };

            let (symbols, energy, steps) = descend(candidate);
            evaluations += steps;
            if !seen.contains(&symbols) {
                offspring = Some((energy, symbols));
                break;
            }
            duplicates += 1;
        }

        let Some(child) = offspring else {
            // every attempt reproduced an ordering already in the population
            break;
        };

        let child_symbols = child.1.clone();
        seen.insert(child_symbols.clone());
        members.push(child);
        members.sort_by(|a, b| a.0.total_cmp(&b.0));
        let dropped = members.pop().expect("population is never empty");
        seen.remove(&dropped.1);
        // the offspring counts as accepted only if it is not the one dropped
        if dropped.1 != child_symbols {
            accepted += 1;
        }

        let best = members[0].0;
        if best < history[history.len() - 1].0 - improvement_tol {
            unsuccessful = 0;
            history.push((best, evaluations));
        } else {
            unsuccessful += 1;
        }
    }

    let energies: Vec<f64> = members.iter().map(|m| m.0).collect();
    let orderings: Vec<Vec<u8>> = members.into_iter().map(|m| m.1).collect();
    Ok((orderings, energies, history, evaluations, generation, accepted, duplicates))
}

#[pymodule]
#[pyo3(name = "_act_rs")]
fn npl_act_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(act_energy, m)?)?;
    m.add_function(wrap_pyfunction!(act_local_optimization, m)?)?;
    m.add_function(wrap_pyfunction!(act_basin_hopping, m)?)?;
    m.add_function(wrap_pyfunction!(act_genetic_algorithm, m)?)?;
    Ok(())
}
