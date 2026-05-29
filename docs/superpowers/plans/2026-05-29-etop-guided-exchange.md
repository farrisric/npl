# ETOP Guided Exchange Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `model="ETOP"` guided-exchange operator so the existing `local_optimization` and `run_basin_hopping` drivers can optimize chemical ordering for an arbitrary number of metals using a trained linear (BayesianRR/TOP) ETOP energy model.

**Architecture:** A new `ETOPExchangeOperator` keeps one `SortedKeyList` per ordered species pair `(i→j)`, ranking atoms by a per-atom flip-gain `g_{i→j}(A) = coef·(f_A^(j) − f_A^(i))` derived from the linear coefficients and the ETOP per-atom feature. The greedy move picks the best i-atom + best j-atom over all pairs; exact accept/reject stays in the drivers, which recompute energy via ETOP's incremental feature update. Driver changes are surgical: a branch in `setup_local_optimization`, a feature-update/rollback branch keyed on `local_env_calculator is None`, and a uniform `guided_exchange` return that removes the binary-only `symbol1_exchange_energies` reach-in.

**Tech Stack:** Python 3.10+, numpy, `sortedcontainers.SortedKeyList`, pytest. Energy model: `BayesianRRCalculator`/`TOPCalculator` (linear, `get_coefficients()`); descriptor: `ExtendedTopologicalFeaturesClassifier` (key `'ETOP'`).

**Spec:** `docs/superpowers/specs/2026-05-29-etop-guided-exchange-design.md`

---

## File structure

- **Create** `npl/optimization/local_optimization/etop_exchange_operator.py` — the `ETOPExchangeOperator`. Sole responsibility: maintain per-species-pair ranked flip-gains and propose guided/hop swaps.
- **Modify** `npl/descriptors/global_feature_classifier.py` — add `ExtendedTopologicalFeaturesClassifier.compute_atom_feature_for_symbol(...)` (pure helper; existing methods untouched).
- **Modify** `npl/optimization/local_optimization/guided_exchange_operator.py` — `guided_exchange` returns the flip estimate as a third value.
- **Modify** `npl/optimization/local_optimization/act_exchange_operator.py` — same `guided_exchange` arity change.
- **Modify** `npl/optimization/local_optimization/local_optimization.py` — register `'ETOP'`, ETOP branch in `setup_local_optimization`, new `update_atomic_features_etop` helper, feature-update/rollback branch in `local_optimization`, 3-value unpack at the `guided_exchange` call site.
- **Modify** `npl/optimization/basin_hopping/basin_hopping.py` — 3-value unpack + drop the `symbol1_exchange_energies` reach-in, feature-update branch in both inner loops.
- **Modify** `tests/test_optimization.py` — add the ETOP test suite (shared `_build_etop()` fixture).

Note on neighborhoods: ETOP per-atom features depend only on an atom's own symbol plus its immediate neighbors' symbols, so a swap perturbs exactly the **one-ring** set `{i1, i2} ∪ neighbors(i1) ∪ neighbors(i2)` (as in `monte_carlo_etop.features_to_update`). The TEC path's two-ring expansion is not needed for ETOP.

---

### Task 1: `compute_atom_feature_for_symbol` on the ETOP classifier

**Files:**
- Modify: `npl/descriptors/global_feature_classifier.py` (add a method to `ExtendedTopologicalFeaturesClassifier`, after `compute_atom_feature`, ~line 236)
- Test: `tests/test_descriptors.py`

- [ ] **Step 1: Write the failing test**

Add to `tests/test_descriptors.py`:

```python
import numpy as np
from npl.core import Nanoparticle
from npl.descriptors import ExtendedTopologicalFeaturesClassifier


def test_compute_atom_feature_for_symbol_matches_current_symbol():
    fc = ExtendedTopologicalFeaturesClassifier(["Au", "Cu", "Pt"])
    p = Nanoparticle()
    p.truncated_octahedron(7, 2, {"Pt": 101, "Au": 50, "Cu": 50})
    fc.compute_feature_vector(p)  # populates sublayer_indices
    idx = p.get_indices_by_symbol("Pt")[0]
    same = fc.compute_atom_feature_for_symbol(p, idx, "Pt")
    assert np.array_equal(same, fc.compute_atom_feature(p, idx))
    # asking for a different symbol changes the feature
    other = fc.compute_atom_feature_for_symbol(p, idx, "Au")
    assert not np.array_equal(same, other)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_descriptors.py::test_compute_atom_feature_for_symbol_matches_current_symbol -v`
Expected: FAIL — `AttributeError: 'ExtendedTopologicalFeaturesClassifier' object has no attribute 'compute_atom_feature_for_symbol'`

- [ ] **Step 3: Add the method**

In `npl/descriptors/global_feature_classifier.py`, inside `ExtendedTopologicalFeaturesClassifier`, immediately after `compute_atom_feature` (before `compute_atom_features`):

```python
    def compute_atom_feature_for_symbol(self, particle, index, symbol):
        """ETOP per-atom feature the atom at ``index`` would have as ``symbol``.

        Geometry (coordination number, sublayer index) and neighbor symbols are
        held fixed; only the atom's own species is hypothetically changed. Relies
        on ``self.sublayer_indices`` being populated (call ``compute_feature_vector``
        once first).
        """
        atom_feature = np.zeros(self.n_features)
        atom_feature[self.layer_types[symbol]] = self.sublayer_indices[index]
        coordination = particle.get_coordination_number(index)
        cn_index = (self.n_bond_features + self.n_sub_layers + coordination +
                    13 * self.symbols.index(symbol))
        for neigh_index in particle.neighbor_list[index]:
            element2 = particle.get_symbol(neigh_index)
            bond_type = tuple(sorted([symbol, element2]))
            atom_feature[self.bond_types[bond_type]] += 0.5
        atom_feature[cn_index] += 1
        return atom_feature
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_descriptors.py::test_compute_atom_feature_for_symbol_matches_current_symbol -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add npl/descriptors/global_feature_classifier.py tests/test_descriptors.py
git commit -m "feat: add ETOP compute_atom_feature_for_symbol helper"
```

---

### Task 2: `ETOPExchangeOperator` — construction, `bind_particle`, flip-gains

**Files:**
- Create: `npl/optimization/local_optimization/etop_exchange_operator.py`
- Test: `tests/test_optimization.py`

- [ ] **Step 1: Add the shared fixture and the gain-formula test**

Add to `tests/test_optimization.py` (top-level, after the existing imports):

```python
import numpy as np
from npl.calculators import BayesianRRCalculator
from npl.optimization.local_optimization.etop_exchange_operator import ETOPExchangeOperator


def _build_etop():
    symbols = ["Au", "Cu", "Pt"]
    fc = ExtendedTopologicalFeaturesClassifier(symbols)
    particle = Nanoparticle()
    particle.truncated_octahedron(7, 2, {"Pt": 101, "Au": 50, "Cu": 50})
    fc.compute_feature_vector(particle)
    coeffs = np.random.default_rng(0).standard_normal(fc.n_features)
    calc = BayesianRRCalculator("ETOP")
    calc.set_coefficients(coeffs)
    return particle, calc, fc


def test_etop_flip_gain_matches_formula():
    particle, calc, fc = _build_etop()
    coeffs = calc.get_coefficients()
    op = ETOPExchangeOperator(coeffs, fc)
    op.bind_particle(particle)

    A = int(particle.get_indices_by_symbol("Pt")[0])
    f_i = fc.compute_atom_feature_for_symbol(particle, A, "Pt")
    f_j = fc.compute_atom_feature_for_symbol(particle, A, "Au")
    expected = float(np.dot(coeffs, f_j) - np.dot(coeffs, f_i))
    assert op.flip_gain[("Pt", "Au")][A] == pytest.approx(expected)
    assert A in op.indices[("Pt", "Au")]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_optimization.py::test_etop_flip_gain_matches_formula -v`
Expected: FAIL — `ModuleNotFoundError`/`ImportError` for `etop_exchange_operator`.

- [ ] **Step 3: Create the operator with `__init__` and `bind_particle`**

Create `npl/optimization/local_optimization/etop_exchange_operator.py`:

```python
import itertools

import numpy as np
from sortedcontainers import SortedKeyList


class ETOPExchangeOperator:
    """Guided chemical-ordering exchange for an arbitrary number of metals.

    For each ordered species pair ``(i, j)`` (i != j) keeps a ``SortedKeyList`` of
    atoms currently of species ``i`` ranked by the flip-gain estimate
    ``g_{i->j}(A) = coef . (f_A^(j) - f_A^(i))``, where ``f_A^(s)`` is A's ETOP
    per-atom feature as if A were species ``s``. A canonical swap of an i-atom A
    with a j-atom B has estimated energy change ``g_{i->j}(A) + g_{j->i}(B)``;
    this is a proposal heuristic (it ignores the i-j cross-bond and shared-neighbor
    coupling). Exact accept/reject is the driver's responsibility.

    Requires a linear energy model: per-atom energy = ``coef . atom_feature``.
    """

    def __init__(self, coefficients, feature_classifier):
        self.coeffs = np.asarray(coefficients)
        self.classifier = feature_classifier
        self.symbols = list(feature_classifier.symbols)
        self.feature_key = feature_classifier.get_feature_key()
        assert len(self.coeffs) == feature_classifier.n_features, (
            "coefficient vector length {} != classifier.n_features {}".format(
                len(self.coeffs), feature_classifier.n_features))

        self.flip_gain = {}     # (sym_i, sym_j) -> {atom_index: gain}
        self.indices = {}       # (sym_i, sym_j) -> SortedKeyList(atom_index)
        self.atom_symbol = {}   # atom_index -> current symbol

    def _pairs(self):
        for i in self.symbols:
            for j in self.symbols:
                if i != j:
                    yield (i, j)

    def _set_gains(self, particle, index):
        sym_i = particle.get_symbol(index)
        self.atom_symbol[index] = sym_i
        base = np.dot(self.coeffs,
                      self.classifier.compute_atom_feature_for_symbol(particle, index, sym_i))
        for sym_j in self.symbols:
            if sym_j == sym_i:
                continue
            f_j = self.classifier.compute_atom_feature_for_symbol(particle, index, sym_j)
            self.flip_gain[(sym_i, sym_j)][index] = float(np.dot(self.coeffs, f_j) - base)
            self.indices[(sym_i, sym_j)].add(index)

    def bind_particle(self, particle):
        self.flip_gain = {pair: {} for pair in self._pairs()}
        self.indices = {
            pair: SortedKeyList(key=lambda idx, p=pair: self.flip_gain[p][idx])
            for pair in self._pairs()
        }
        self.atom_symbol = {}
        for index in particle.get_indices():
            self._set_gains(particle, index)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_optimization.py::test_etop_flip_gain_matches_formula -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add npl/optimization/local_optimization/etop_exchange_operator.py tests/test_optimization.py
git commit -m "feat: ETOPExchangeOperator construction and flip-gain binding"
```

---

### Task 3: `guided_exchange`

**Files:**
- Modify: `npl/optimization/local_optimization/etop_exchange_operator.py`
- Test: `tests/test_optimization.py`

- [ ] **Step 1: Write the failing test**

Add to `tests/test_optimization.py`:

```python
def test_etop_guided_exchange_swaps_unlike_species():
    particle, calc, fc = _build_etop()
    op = ETOPExchangeOperator(calc.get_coefficients(), fc)
    op.bind_particle(particle)
    sym_a = particle.get_symbol  # alias for readability

    A, B, flip = op.guided_exchange(particle)
    # the two atoms were of different species before the swap
    assert op.atom_symbol[A] != op.atom_symbol[B]
    # after the in-place swap they carry each other's old species
    assert sym_a(A) == op.atom_symbol[B]
    assert sym_a(B) == op.atom_symbol[A]
    assert isinstance(flip, float)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_optimization.py::test_etop_guided_exchange_swaps_unlike_species -v`
Expected: FAIL — `AttributeError: 'ETOPExchangeOperator' object has no attribute 'guided_exchange'`

- [ ] **Step 3: Implement `guided_exchange`**

Add to `ETOPExchangeOperator` (after `bind_particle`):

```python
    def guided_exchange(self, particle):
        best = None  # (gain, A, B)
        for i, j in itertools.combinations(self.symbols, 2):
            la, lb = self.indices[(i, j)], self.indices[(j, i)]
            if la and lb:
                a, b = la[0], lb[0]
                gain = self.flip_gain[(i, j)][a] + self.flip_gain[(j, i)][b]
                if best is None or gain < best[0]:
                    best = (gain, a, b)
        if best is None:
            raise ValueError(
                "ETOP guided exchange needs at least two species present in the particle.")
        gain, index_a, index_b = best
        particle.swap_symbols([(index_a, index_b)])
        return index_a, index_b, gain
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_optimization.py::test_etop_guided_exchange_swaps_unlike_species -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add npl/optimization/local_optimization/etop_exchange_operator.py tests/test_optimization.py
git commit -m "feat: ETOPExchangeOperator.guided_exchange"
```

---

### Task 4: `update` (incremental refresh) — the key invariant

**Files:**
- Modify: `npl/optimization/local_optimization/etop_exchange_operator.py`
- Test: `tests/test_optimization.py`

- [ ] **Step 1: Write the failing test**

Add to `tests/test_optimization.py`:

```python
def test_etop_update_matches_fresh_bind():
    particle, calc, fc = _build_etop()
    coeffs = calc.get_coefficients()
    op = ETOPExchangeOperator(coeffs, fc)
    op.bind_particle(particle)

    for _ in range(5):
        a, b, _flip = op.guided_exchange(particle)
        neighborhood = {a, b}
        neighborhood |= set(particle.neighbor_list[a])
        neighborhood |= set(particle.neighbor_list[b])
        fc.update_feature_vector(particle, neighborhood)
        op.update(particle, neighborhood, [a, b])

    fresh = ETOPExchangeOperator(coeffs, fc)
    fresh.bind_particle(particle)
    for pair in op.flip_gain:
        assert op.flip_gain[pair].keys() == fresh.flip_gain[pair].keys()
        for idx in op.flip_gain[pair]:
            assert op.flip_gain[pair][idx] == pytest.approx(fresh.flip_gain[pair][idx])
        assert list(op.indices[pair]) == list(fresh.indices[pair])
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_optimization.py::test_etop_update_matches_fresh_bind -v`
Expected: FAIL — `AttributeError: 'ETOPExchangeOperator' object has no attribute 'update'`

- [ ] **Step 3: Implement `update`**

Add to `ETOPExchangeOperator`. The remove phase MUST run before recomputation: `SortedKeyList.remove` re-reads the key, so entries must be removed while their stored gain is still the old value.

```python
    def update(self, particle, indices, exchange_indices):
        # exchange_indices is accepted for interface parity with the binary
        # operators; ETOP refreshes every atom in `indices` from self.atom_symbol.
        for index in indices:
            sym_old = self.atom_symbol[index]
            for sym_j in self.symbols:
                if sym_j != sym_old:
                    self.indices[(sym_old, sym_j)].remove(index)
                    del self.flip_gain[(sym_old, sym_j)][index]
        for index in indices:
            self._set_gains(particle, index)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_optimization.py::test_etop_update_matches_fresh_bind -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add npl/optimization/local_optimization/etop_exchange_operator.py tests/test_optimization.py
git commit -m "feat: ETOPExchangeOperator.update with incremental==fresh invariant"
```

---

### Task 5: `basin_hop_step`

**Files:**
- Modify: `npl/optimization/local_optimization/etop_exchange_operator.py`
- Test: `tests/test_optimization.py`

- [ ] **Step 1: Write the failing test**

Add to `tests/test_optimization.py`:

```python
def test_etop_basin_hop_step_swaps_unlike_species():
    particle, calc, fc = _build_etop()
    op = ETOPExchangeOperator(calc.get_coefficients(), fc)
    op.bind_particle(particle)
    a, b = op.basin_hop_step(particle)
    assert op.atom_symbol[a] != op.atom_symbol[b]
    assert particle.get_symbol(a) == op.atom_symbol[b]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_optimization.py::test_etop_basin_hop_step_swaps_unlike_species -v`
Expected: FAIL — `AttributeError: ... 'basin_hop_step'`

- [ ] **Step 3: Implement `basin_hop_step`**

A hop should perturb away from the greedy minimum: prefer the least-uphill candidate (smallest positive estimated ΔE); if every candidate is still downhill, take the least beneficial (largest gain). Add to `ETOPExchangeOperator`:

```python
    def basin_hop_step(self, particle):
        candidates = []  # (gain, A, B)
        for i, j in itertools.combinations(self.symbols, 2):
            la, lb = self.indices[(i, j)], self.indices[(j, i)]
            if la and lb:
                a, b = la[0], lb[0]
                candidates.append((self.flip_gain[(i, j)][a] + self.flip_gain[(j, i)][b], a, b))
        if not candidates:
            raise ValueError(
                "ETOP basin hop needs at least two species present in the particle.")
        uphill = [c for c in candidates if c[0] > 0]
        gain, index_a, index_b = (min(uphill, key=lambda c: c[0]) if uphill
                                  else max(candidates, key=lambda c: c[0]))
        particle.swap_symbols([(index_a, index_b)])
        return index_a, index_b
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_optimization.py::test_etop_basin_hop_step_swaps_unlike_species -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add npl/optimization/local_optimization/etop_exchange_operator.py tests/test_optimization.py
git commit -m "feat: ETOPExchangeOperator.basin_hop_step"
```

---

### Task 6: Uniform `guided_exchange` return for the binary operators + call sites

This removes the binary-only `symbol1_exchange_energies` reach-in from `basin_hopping`, so all operators share one return contract `(index1, index2, flip)`. TOP/ACT behavior is unchanged.

**Files:**
- Modify: `npl/optimization/local_optimization/guided_exchange_operator.py:24-29`
- Modify: `npl/optimization/local_optimization/act_exchange_operator.py:32-37`
- Modify: `npl/optimization/local_optimization/local_optimization.py:65`
- Modify: `npl/optimization/basin_hopping/basin_hopping.py:36-39`
- Test: existing `tests/test_optimization.py` TOP/ACT tests

- [ ] **Step 1: Confirm the guard tests exist and pass now**

Run: `pytest tests/test_optimization.py -k "TOP or ACT" -v`
Expected: PASS (these are the regression guard for this task).

- [ ] **Step 2: Update `GuidedExchangeOperator.guided_exchange`**

Replace the body of `guided_exchange` in `guided_exchange_operator.py`:

```python
    def guided_exchange(self, particle):
        symbol1_index = self.symbol1_indices[0]
        symbol2_index = self.symbol2_indices[0]

        flip = (self.symbol1_exchange_energies[symbol1_index]
                + self.symbol2_exchange_energies[symbol2_index])
        particle.swap_symbols([(symbol1_index, symbol2_index)])
        return symbol1_index, symbol2_index, flip
```

- [ ] **Step 3: Update `ACTExchangeOperator.guided_exchange`**

Replace the body of `guided_exchange` in `act_exchange_operator.py`:

```python
    def guided_exchange(self, particle):
        symbol1_index = self.symbol1_indices[0]
        symbol2_index = self.symbol2_indices[0]

        flip = (self.symbol1_exchange_energies[symbol1_index]
                + self.symbol2_exchange_energies[symbol2_index])
        particle.swap_symbols([(symbol1_index, symbol2_index)])
        return symbol1_index, symbol2_index, flip
```

- [ ] **Step 4: Update the `local_optimization` call site**

In `local_optimization.py`, line 65, change:

```python
        index1, index2 = exchange_operator.guided_exchange(start_particle)
```
to:
```python
        index1, index2, _flip = exchange_operator.guided_exchange(start_particle)
```

- [ ] **Step 5: Update the `basin_hopping` call site (drop the reach-in)**

In `basin_hopping.py`, replace lines 36–39:

```python
            index1, index2 = exchange_operator.guided_exchange(start_particle)

            flip = (exchange_operator.symbol1_exchange_energies[index1]
                    + exchange_operator.symbol2_exchange_energies[index2])
```
with:
```python
            index1, index2, flip = exchange_operator.guided_exchange(start_particle)
```

- [ ] **Step 6: Run the guard tests**

Run: `pytest tests/test_optimization.py -k "TOP or ACT" -v`
Expected: PASS (4 tests).

- [ ] **Step 7: Commit**

```bash
git add npl/optimization/local_optimization/guided_exchange_operator.py \
        npl/optimization/local_optimization/act_exchange_operator.py \
        npl/optimization/local_optimization/local_optimization.py \
        npl/optimization/basin_hopping/basin_hopping.py
git commit -m "refactor: guided_exchange returns flip estimate uniformly"
```

---

### Task 7: ETOP wiring in `local_optimization` + error path

**Files:**
- Modify: `npl/optimization/local_optimization/local_optimization.py`
- Test: `tests/test_optimization.py`

- [ ] **Step 1: Write the failing tests**

Add to `tests/test_optimization.py`:

```python
def test_etop_local_optimization_descends_and_preserves_composition():
    particle, calc, fc = _build_etop()
    before = dict(particle.get_stoichiometry())
    result, accepted = local_optimization(particle, calc, None, model="ETOP")
    energies = [e for e, _ in accepted]
    assert all(energies[k + 1] <= energies[k] + 1e-9 for k in range(len(energies) - 1))
    assert accepted[-1][0] <= accepted[0][0] + 1e-6
    assert dict(result.get_stoichiometry()) == before


def test_etop_requires_linear_calculator():
    particle, _calc, _fc = _build_etop()

    class _NoCoeffCalc:
        def get_energy_key(self):
            return "X"

        def compute_energy(self, p):
            p.set_energy("X", 0.0)

    with pytest.raises(ValueError):
        local_optimization(particle, _NoCoeffCalc(), None, model="ETOP")
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_optimization.py -k etop_local_optimization -v`
Expected: FAIL — `KeyError: 'ETOP'` (model not registered / no ETOP branch).

- [ ] **Step 3: Add imports and register the model**

In `local_optimization.py`, add to the imports at the top:

```python
from npl.descriptors.global_feature_classifier import ExtendedTopologicalFeaturesClassifier
from npl.optimization.local_optimization.etop_exchange_operator import ETOPExchangeOperator
```

And extend the registry:

```python
EXCHANGE_OPERATORS = {'TOP': GuidedExchangeOperator, 'ACT': ACTExchangeOperator,
                      'ETOP': ETOPExchangeOperator}
```

- [ ] **Step 4: Add the ETOP branch to `setup_local_optimization`**

Insert at the very start of the function body (before `local_env_calculator = ...`):

```python
    symbols = start_particle.get_all_symbols()

    if model == 'ETOP':
        if not hasattr(energy_calculator, 'get_coefficients'):
            raise ValueError(
                "ETOP guided exchange requires a linear energy model exposing "
                "get_coefficients() (e.g. BayesianRRCalculator or TOPCalculator).")
        if local_feature_classifier is None:
            local_feature_classifier = ExtendedTopologicalFeaturesClassifier(symbols)
        local_feature_classifier.compute_feature_vector(start_particle)
        energy_calculator.compute_energy(start_particle)
        exchange_operator = ETOPExchangeOperator(
            energy_calculator.get_coefficients(), local_feature_classifier)
        exchange_operator.bind_particle(start_particle)
        energy_key = energy_calculator.get_energy_key()
        return energy_key, None, local_feature_classifier, exchange_operator
```

Then delete the now-duplicate `symbols = start_particle.get_all_symbols()` line that previously opened the TOP/ACT path (it has moved above the branch).

- [ ] **Step 5: Add the `update_atomic_features_etop` helper**

In `local_optimization.py`, after the existing `update_atomic_features` function:

```python
def update_atomic_features_etop(index1, index2, feature_classifier, particle):
    neighborhood = {index1, index2}
    neighborhood = neighborhood.union(particle.neighbor_list[index1])
    neighborhood = neighborhood.union(particle.neighbor_list[index2])

    old_atom_features, change = feature_classifier.update_feature_vector(particle, neighborhood)

    return particle, neighborhood, old_atom_features, change
```

- [ ] **Step 6: Branch the feature-update and rollback in `local_optimization`**

Replace the `while True` loop body (current lines 63–89) with:

```python
    while True:
        step += 2
        index1, index2, _flip = exchange_operator.guided_exchange(start_particle)
        exchanged_indices = [index1, index2]

        if local_env_calculator is None:
            start_particle, neighborhood, old_atom_features, change = \
                update_atomic_features_etop(index1, index2, local_feature_classifier,
                                            start_particle)
        else:
            start_particle, neighborhood = update_atomic_features(
                index1, index2, local_env_calculator, local_feature_classifier, start_particle)
        exchange_operator.update(start_particle, neighborhood, exchanged_indices)

        energy_calculator.compute_energy(start_particle)
        new_energy = start_particle.get_energy(energy_key)

        if new_energy < start_energy:
            start_energy = new_energy
            accepted_energies.append((new_energy, step))
        else:
            accepted_energies.append((start_energy, step))

            # roll back last exchange and make sure features and environments are up-to-date
            start_particle.swap_symbols([(index1, index2)])
            start_particle.set_energy(energy_key, start_energy)
            if local_env_calculator is None:
                local_feature_classifier.downgrade_feature_vector(
                    start_particle, neighborhood, old_atom_features, change)
            else:
                for index in neighborhood:
                    local_env_calculator.compute_local_environment(start_particle, index)
                    local_feature_classifier.compute_atom_feature(start_particle, index)

            break
```

- [ ] **Step 7: Run tests to verify they pass**

Run: `pytest tests/test_optimization.py -k "etop_local_optimization or etop_requires" -v`
Expected: PASS (2 tests).

- [ ] **Step 8: Commit**

```bash
git add npl/optimization/local_optimization/local_optimization.py tests/test_optimization.py
git commit -m "feat: wire model=ETOP into local_optimization"
```

---

### Task 8: ETOP wiring in `run_basin_hopping`

**Files:**
- Modify: `npl/optimization/basin_hopping/basin_hopping.py`
- Test: `tests/test_optimization.py`

- [ ] **Step 1: Write the failing test**

Add to `tests/test_optimization.py`:

```python
def test_etop_basin_hopping_runs():
    particle, calc, fc = _build_etop()
    best, lowest_energies, flip_list = run_basin_hopping(
        particle, calc, None, n_hopping_attempts=2, n_hops=2, model="ETOP")
    assert best is not None
    assert len(lowest_energies) >= 1
    assert lowest_energies[-1][0] <= lowest_energies[0][0] + 1e-6
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_optimization.py::test_etop_basin_hopping_runs -v`
Expected: FAIL — `AttributeError` (`update_atomic_features` calls `local_env_calculator.compute_local_environment` but `local_env_calculator` is `None` for ETOP).

- [ ] **Step 3: Import the ETOP helper**

In `basin_hopping.py`, change the import line:

```python
from npl.optimization.local_optimization import update_atomic_features
```
to:
```python
from npl.optimization.local_optimization import update_atomic_features
from npl.optimization.local_optimization import update_atomic_features_etop
```

- [ ] **Step 4: Branch the feature update in the local-opt inner loop**

In `basin_hopping.py`, replace the block at lines 43–46:

```python
            start_particle, neighborhood = update_atomic_features(index1, index2,
                                                                  local_env_calculator,
                                                                  local_feature_classifier,
                                                                  start_particle)
```
with:
```python
            if local_env_calculator is None:
                start_particle, neighborhood, _old, _chg = update_atomic_features_etop(
                    index1, index2, local_feature_classifier, start_particle)
            else:
                start_particle, neighborhood = update_atomic_features(
                    index1, index2, local_env_calculator, local_feature_classifier,
                    start_particle)
```

- [ ] **Step 5: Branch the feature update in the hop loop**

In `basin_hopping.py`, replace the block at lines 82–85:

```python
            start_particle, neighborhood = update_atomic_features(index1, index2,
                                                                  local_env_calculator,
                                                                  local_feature_classifier,
                                                                  start_particle)
```
with:
```python
            if local_env_calculator is None:
                start_particle, neighborhood, _old, _chg = update_atomic_features_etop(
                    index1, index2, local_feature_classifier, start_particle)
            else:
                start_particle, neighborhood = update_atomic_features(
                    index1, index2, local_env_calculator, local_feature_classifier,
                    start_particle)
```

- [ ] **Step 6: Run test to verify it passes**

Run: `pytest tests/test_optimization.py::test_etop_basin_hopping_runs -v`
Expected: PASS

- [ ] **Step 7: Commit**

```bash
git add npl/optimization/basin_hopping/basin_hopping.py tests/test_optimization.py
git commit -m "feat: wire model=ETOP into run_basin_hopping"
```

---

### Task 9: Full suite + lint

**Files:** none (verification only)

- [ ] **Step 1: Run the whole optimization + descriptor suite**

Run: `pytest tests/test_optimization.py tests/test_descriptors.py -v`
Expected: PASS (all TOP/ACT tests plus the new ETOP tests).

- [ ] **Step 2: Run the entire test suite**

Run: `pytest`
Expected: PASS (no regressions).

- [ ] **Step 3: Lint**

Run: `flake8`
Expected: no errors (max line 100; mind line length in the new operator and the branch blocks).

- [ ] **Step 4: Commit any lint fixes**

```bash
git add -A
git commit -m "style: flake8 fixes for ETOP guided exchange"
```

---

## Self-review notes

- **Spec coverage:** operator core (Tasks 2–5), `model="ETOP"` drop-in for both drivers (Tasks 7–8), error handling — non-linear calculator (Task 7, `test_etop_requires_linear_calculator`), coefficient-length assertion (Task 2 `__init__`), single-species guard (Tasks 3/5 `ValueError`), `environment_energies` ignored for ETOP (passed as `None` in all ETOP tests). Tests map to spec §Testing: gain-formula correctness (Task 2), incremental==fresh (Task 4), monotonic descent + composition (Task 7), basin-hopping smoke (Task 8), error path (Task 7).
- **Spec correction:** spec test #1 originally claimed `flip_estimate == exact ΔE` for non-adjacent atoms; that is false because the gain weights each bond at 0.5 (own half only). Replaced with the white-box gain-formula test in Task 2 and corrected in the spec doc.
- **Type/name consistency:** operator exposes `flip_gain`, `indices`, `atom_symbol`, `guided_exchange -> (a, b, flip)`, `basin_hop_step -> (a, b)`, `update(particle, indices, exchange_indices)`, `bind_particle(particle)`; helper `update_atomic_features_etop(index1, index2, feature_classifier, particle) -> (particle, neighborhood, old_atom_features, change)`; classifier method `compute_atom_feature_for_symbol(particle, index, symbol)`. These names are used identically across tasks.
