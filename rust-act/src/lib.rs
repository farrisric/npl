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

#[pymodule]
fn npl_act_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(act_energy, m)?)?;
    m.add_function(wrap_pyfunction!(act_local_optimization, m)?)?;
    m.add_function(wrap_pyfunction!(act_basin_hopping, m)?)?;
    Ok(())
}
