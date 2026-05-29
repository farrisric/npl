# ETOP guided exchange for k ≥ 2 metals (`model="ETOP"`)

**Date:** 2026-05-29
**Status:** Design approved, pending implementation plan

## Problem

`npl.optimization` provides *guided* chemical-ordering optimization (basin hopping +
greedy local optimization) only for **two** elements. The exchange operators
(`GuidedExchangeOperator` for `model="TOP"`, `ACTExchangeOperator` for `model="ACT"`)
encode a binary flip: each atom carries a single scalar `env_energy_difference`, and the
best move is the top of two sorted lists (`symbol1_indices`, `symbol2_indices`). This
structure cannot represent ≥3 species, so multi-metal ordering currently has no guided
path and must use the Metropolis MC in `npl.monte_carlo`.

We want a drop-in `model="ETOP"` exchange operator that works in the **existing**
`local_optimization` and `run_basin_hopping` drivers for an arbitrary number of metals,
driven by a **trained linear (BayesianRR) ETOP energy model**.

## Why a linear ETOP model makes this exact

`BayesianRRCalculator` is `fit_intercept=False` and computes
`E = np.dot(coef_, feature_vector)`. The ETOP global feature vector is the sum of per-atom
features (`ExtendedTopologicalFeaturesClassifier.compute_feature_vector` →
`atom_features.sum(axis=0)`). Therefore the energy decomposes **exactly** per atom:

```
E_total = coef_ · Σ_A f_A = Σ_A (coef_ · f_A)
```

The guided heuristic and the exact energy share one coefficient vector
(`energy_calculator.get_coefficients()` → `ridge.coef_`). No intercept, no separate
environment-energy table is needed.

## Scope

- **In:** new `ETOPExchangeOperator`; registration as `model="ETOP"`; the surgical driver
  changes needed for it to run inside `local_optimization` and `run_basin_hopping`; tests.
- **Out:** non-linear energy models (GPR/kernel) — rejected at setup with a clear error;
  changes to the MC path; changes to TOP/ACT numerical behavior.

## Core algorithm

`ETOPExchangeOperator` keeps, for each **ordered species pair** `(i → j)` with `i ≠ j`, a
`SortedKeyList` of atoms currently of species `i`, keyed by a flip-gain estimate:

```
g_{i→j}(A) = coef_ · ( f_A^(j) − f_A^(i) )
```

where `f_A^(s)` is A's ETOP atom-feature computed as if A were species `s`, holding A's
coordination number, neighbor symbols, and sublayer index fixed. This depends only on A's
local neighborhood, so it is cheap. With `k` species there are `k(k−1)` directed lists.

A canonical swap exchanges an i-atom A with a j-atom B (A becomes j, B becomes i). Its
estimated energy change is

```
ΔE_est = g_{i→j}(A) + g_{j→i}(B)
```

For each unordered pair `{i, j}` the best candidate is `top(i→j) + top(j→i)`; the proposed
move is the globally most-negative candidate over all pairs. This is a **proposal
heuristic**: it ignores the i–j cross-bond (if A, B are neighbors) and shared-neighbor
coupling — exactly as the binary operator does. Correctness comes from the driver, which
recomputes the exact energy after the swap and rolls back if it did not improve (local
optimization) or accepts per the basin-hopping logic.

### Operator interface

`etop_exchange_operator.py`:

- `__init__(self, coefficients, feature_classifier)` — stores `coef_`,
  `symbols = feature_classifier.symbols` (sorted), `feature_key = 'ETOP'`, and a per-atom
  current-symbol map. Validates `len(coefficients) == feature_classifier.n_features`.
- `bind_particle(self, particle)` — builds all `k(k−1)` directed sorted lists.
- `guided_exchange(self, particle)` — returns `(index_a, index_b, flip_estimate)`;
  performs the swap via `particle.swap_symbols`. Skips directed lists that are empty
  (species absent). Raises if no swappable pair exists (single-element particle).
- `basin_hop_step(self, particle)` — proposes a non-greedy, least-uphill move to escape a
  local minimum, mirroring the scan structure of the binary operator's `basin_hop_step`.
  Returns `(index_a, index_b)`.
- `update(self, particle, exchanged_indices, neighborhood)` — for each atom in the
  recomputed `neighborhood`: remove it from the directed lists keyed by its previous
  symbol, recompute its flip-gains for the new symbol, re-add. Mirrors the binary
  operator's `update` remove/recompute/add pattern.

### Helper on the classifier

Add `ExtendedTopologicalFeaturesClassifier.compute_atom_feature_for_symbol(particle,
index, symbol)` returning the atom-feature vector A would have as `symbol`. This is the
existing `compute_atom_feature` logic with the element taken from the argument instead of
`particle.get_symbol(index)`. Additive; existing behavior unchanged. Relies on
`self.sublayer_indices` already being populated (it is, after `compute_feature_vector`).

## Driver integration (surgical, Approach 1)

`local_optimization/local_optimization.py`:

- `EXCHANGE_OPERATORS['ETOP'] = ETOPExchangeOperator`.
- `setup_local_optimization`: add an `if model == 'ETOP'` branch that
  - uses `local_feature_classifier or ExtendedTopologicalFeaturesClassifier(symbols)`,
  - computes the feature vector and energy,
  - reads `coefficients = energy_calculator.get_coefficients()`, raising
    `ValueError("ETOP guided exchange requires a linear energy model exposing "
    "get_coefficients()")` if the calculator has no such method,
  - builds and binds `ETOPExchangeOperator(coefficients, classifier)`,
  - returns `local_env_calculator = None`.
  The existing TOP/ACT branch is untouched. The `environment_energies` argument is unused
  (and documented as ignored) for ETOP.
- Feature update / rollback: branch on `local_env_calculator is None`. For ETOP, use the
  classifier's incremental `update_feature_vector` over a one-ring neighborhood
  (`{i1, i2} ∪ neighbors(i1) ∪ neighbors(i2)`, as in `monte_carlo_etop.features_to_update`)
  and `downgrade_feature_vector` for rollback. The TEC path (two-ring, via
  `local_env_calculator`) is unchanged.

`guided_exchange` return uniformity:

- Both `GuidedExchangeOperator.guided_exchange` and `ETOPExchangeOperator.guided_exchange`
  return `(index1, index2, flip_estimate)`. For the binary operator, `flip_estimate =
  symbol1_exchange_energies[i] + symbol2_exchange_energies[j]` (the value
  `basin_hopping.py` computes today at lines 42–43).
- Call sites updated: `local_optimization` line 65 (`index1, index2, _ = ...`) and
  `basin_hopping` lines 40–43 (use the returned `flip`, drop the
  `symbol1_exchange_energies[...]` reach-in). `flip_energy_list` content is preserved.

`basin_hopping/basin_hopping.py`:

- Feature-update and rollback call the same `local_env_calculator is None` branch.
- Final best-particle rebuild (lines 99–105): for ETOP (`local_env_calculator is None`)
  skip `compute_local_environments` and just call `compute_feature_vector(best_particle)`,
  which recomputes ETOP features from symbols.

## Error handling

- `model='ETOP'` with a calculator lacking `get_coefficients` → `ValueError` (message
  above).
- `len(coefficients) != classifier.n_features` → assertion with a descriptive message.
- Single-element particle (no swappable species pair) → clear error from
  `guided_exchange` / `bind_particle`.
- `environment_energies` ignored for ETOP — documented, not an error.

## Testing (`tests/test_optimization.py`)

1. **Linear-decomposition exactness.** Construct a particle and a linear ETOP calculator;
   choose a swap where A and B are non-adjacent and share no neighbors. Assert the
   operator's `flip_estimate` equals the exact ΔE from a full energy recompute (the only
   terms that could differ — cross-bond and shared neighbors — are absent).
2. **Incremental == fresh.** After a sequence of guided swaps, build a fresh operator via
   `bind_particle` on the current particle and assert its directed-list contents and
   flip-gains match the incrementally `update`d operator. This is the primary invariant
   guarding the `update` logic.
3. **Monotonic descent + composition preserved.** `local_optimization(model='ETOP')`
   yields a non-increasing accepted-energy sequence, terminates at a local minimum, and
   leaves per-species counts unchanged.
4. **Basin-hopping smoke.** `run_basin_hopping(model='ETOP')` runs and returns a best
   energy ≤ the start energy.
5. **Error path.** A calculator without `get_coefficients` (e.g. a GPR stub) raises
   `ValueError` at setup.

## Notes on equivalence

There is **no** numeric reduction to the binary `GuidedExchangeOperator` at `k = 2`: the
two share the per-pair-sorted-list *structure* (two directed lists), but they consume
different descriptors (ETOP vs. TEC), so their energies are not comparable. The real
correctness invariants are tests 1 (exactness) and 2 (incremental == fresh).
