# NPL JOSS-Readiness Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make NanoParticleLibrary (`npl`) JOSS-submittable: it imports and runs on a current scientific-Python stack, dead code and the migrated-out Monte Carlo module are gone, an automated pytest suite runs in CI, packaging is modern (`pyproject.toml`), and a JOSS paper is drafted.

**Architecture:** Modernize in verified vertical slices. Six phases, each with an exit criterion that must pass before the next begins: (1) restore importability, (2) remove `npl.monte_carlo` and redraw the boundary, (3) dead-code + latent-bug sweep, (4) tests + CI, (5) packaging, (6) paper + repo essentials.

**Tech Stack:** Python 3.10–3.13, NumPy 2.x, SciPy 1.17, ASE 3.28, scikit-learn 1.8, pytest, GitHub Actions, setuptools/PEP 621.

**Branch:** all work on `joss-readiness` (already created). **Commits:** plain messages, no `Co-Authored-By`/AI-authorship trailer.

**Reference design:** `docs/superpowers/specs/2026-05-29-npl-joss-readiness-design.md`. Sibling repo `mcpy` (`/home/energystorage/projects/mcpy`) is the quality template for `pyproject.toml`, `CONTRIBUTING.md`, and `paper/`.

---

## File Structure

Created:
- `tests/test_imports.py`, `tests/test_core.py`, `tests/test_descriptors.py`, `tests/test_calculators.py`, `tests/test_optimization.py`, `tests/conftest.py`
- `examples/top_energy_evaluation.py` (runnable, CI-tested flagship example)
- `pyproject.toml`
- `paper/paper.md`, `paper/paper.bib`
- `CONTRIBUTING.md`

Modified:
- `npl/descriptors/local_environment_calculator.py` (remove dead `SOAPCalculator` + `sph_harm` import)
- `npl/core/nanoparticle.py` (`np.int` → `int`)
- `npl/calculators/energy_calculator.py` (`set_calculator` → `.calc`; remove dead classes/fn; drop commented code)
- `npl/descriptors/__init__.py`, `npl/optimization/__init__.py` (fix `__all__`)
- `npl/optimization/basin_hopping/basin_hopping.py` (fix `flip` expression bug; thread `model` selector)
- `npl/optimization/local_optimization/local_optimization.py` (add `model` selector for TOP/ACT operator)
- `npl/__init__.py` (version → 1.1.0)
- `.github/workflows/test.yml`, `.github/workflows/publish.yml`
- `README.md`, `docs/*` (de-reference MC; Python version)
- `CITATION.cff`

Deleted:
- `npl/monte_carlo/` (whole tree)
- `setup.py`, `data/params.json`
- `test/test_core.py` (empty; replaced by `tests/`)
- `docs/npl.monte_carlo.rst`, `docs/npl.monte_carlo.ensembles.rst`

Renamed:
- `npl/optimization/genetic_algoritm/` → `npl/optimization/genetic_algorithm/`
- `npl/optimization/local_optimization/garbage_exchange_operator.py` → `act_exchange_operator.py` (class `GuidedExchangeOperator` → `ACTExchangeOperator`)

---

# Phase 1 — Restore importability on the current stack

**Exit criterion:** `python -c "import npl, npl.core, npl.descriptors, npl.calculators, npl.optimization, npl.visualize, npl.utils"` exits 0.

### Task 1.1: Remove the dead `SOAPCalculator` and its `sph_harm` import

`SOAPCalculator` is unused anywhere in the repo (`grep -rn SOAPCalculator npl examples docs` → only its definition) and is the sole reason for the `from scipy.special import sph_harm` import that breaks the whole `npl.descriptors` import chain on SciPy ≥ 1.15. Deleting it both restores importability and removes dead code.

**Files:**
- Modify: `npl/descriptors/local_environment_calculator.py`

- [ ] **Step 1: Delete the `sph_harm` import and the entire `SOAPCalculator` class**

In `npl/descriptors/local_environment_calculator.py`, remove line 3 (`from scipy.special import sph_harm`) and the whole `class SOAPCalculator(LocalEnvironmentCalculator):` block (lines 27–108). Keep `import numpy as np` and `import copy` (still used by `NeighborCountingEnvironmentCalculator`). The file should retain only `LocalEnvironmentCalculator` and `NeighborCountingEnvironmentCalculator`.

- [ ] **Step 2: Verify the module imports**

Run: `python -c "import npl.descriptors.local_environment_calculator as m; print(sorted(n for n in dir(m) if not n.startswith('_')))"`
Expected: prints a list containing `LocalEnvironmentCalculator` and `NeighborCountingEnvironmentCalculator`, no `SOAPCalculator`, no ImportError.

### Task 1.2: Fix removed NumPy alias `np.int`

**Files:**
- Modify: `npl/core/nanoparticle.py:44`

- [ ] **Step 1: Replace `np.int` with builtin `int`**

Change `excess_atoms = excess_atoms.astype(np.int)` to `excess_atoms = excess_atoms.astype(int)`.

- [ ] **Step 2: Verify**

Run: `python -c "import npl.core.nanoparticle"`
Expected: exit 0, no `AttributeError: module 'numpy' has no attribute 'int'`.

### Task 1.3: Fix removed ASE API `set_calculator`

**Files:**
- Modify: `npl/calculators/energy_calculator.py:104`

- [ ] **Step 1: Replace `set_calculator` with the `.calc` attribute**

Change `atoms.set_calculator(EMT())` to `atoms.calc = EMT()`.

- [ ] **Step 2: Verify imports of the calculators package**

Run: `python -c "import npl.calculators"`
Expected: exit 0.

### Task 1.4: Full import sweep + commit

- [ ] **Step 1: Confirm every public submodule imports**

Run:
```bash
python -c "import npl, npl.core, npl.descriptors, npl.calculators, npl.optimization, npl.visualize, npl.utils; print('all import OK')"
```
Expected: `all import OK`. If a new error surfaces, fix the specific removed/renamed API (search the traceback's file/line), then re-run until clean.

- [ ] **Step 2: Commit**

```bash
git add npl/descriptors/local_environment_calculator.py npl/core/nanoparticle.py npl/calculators/energy_calculator.py
git commit -m "fix: restore importability on numpy 2 / scipy 1.17 / ASE 3.28

Remove unused SOAPCalculator (dropped its removed scipy.special.sph_harm
import), replace np.int with int, and atoms.set_calculator with atoms.calc."
```

---

# Phase 2 — Remove `npl.monte_carlo`; redraw the boundary

**Exit criterion:** no references to `npl.monte_carlo` remain in `npl/`, `examples/`, `docs/`, or `README.md`; `examples/top_energy_evaluation.py` runs and prints a finite energy.

### Task 2.1: Delete the `npl.monte_carlo` package

`npl.monte_carlo` has zero internal importers (`grep -rn monte_carlo npl --include=*.py | grep -v '^npl/monte_carlo/'` → empty). It is the MC code that moved to the `mcpy` repo. The README and several notebooks reference it; those are handled in later tasks.

**Files:**
- Delete: `npl/monte_carlo/` (entire directory)

- [ ] **Step 1: Remove the directory**

Run: `git rm -r npl/monte_carlo`

- [ ] **Step 2: Verify package still imports**

Run: `python -c "import npl, npl.optimization, npl.calculators; print('ok')"`
Expected: `ok` (nothing imported `npl.monte_carlo`).

### Task 2.2: Fix broken `__all__` exports

Two packages export names that are never bound: `npl/descriptors/__init__.py` lists `"LayererTopologicalDescriptors"` (defined in `global_feature_classifier.py` but never imported into the package) and `"LocalEnvironmentCalculator"` (this one *is* imported, keep it); `npl/optimization/__init__.py` has a missing comma producing `"local_optimization" "run_basin_hopping"` (one concatenated string) and lists `GuidedSearch`/`MCSearch`/`GASearch` which are real.

**Files:**
- Modify: `npl/descriptors/__init__.py`
- Modify: `npl/optimization/__init__.py`

- [ ] **Step 1: Make `descriptors/__all__` match what is actually exported**

Replace the `__all__` block at the bottom of `npl/descriptors/__init__.py` with one that lists only names imported in that file:

```python
__all__ = [
    "GlobalFeatureClassifier",
    "TopologicalFeatureClassifier",
    "ExtendedTopologicalFeaturesClassifier",
    "AtomicCoordinationTypes",
    "CoordinationFeatureClassifier",
    "LocalEnvironmentFeatureClassifier",
    "CoordinationNumberClassifier",
    "TopologicalEnvironmentClassifier",
    "LocalEnvironmentCalculator",
    "NeighborCountingEnvironmentCalculator",
]
```

- [ ] **Step 2: Fix the comma bug in `optimization/__all__`**

Replace the `__all__` block in `npl/optimization/__init__.py` with:

```python
__all__ = [
    "GOSearch",
    "MCSearch",
    "GASearch",
    "GuidedSearch",
    "local_optimization",
    "run_basin_hopping",
]
```

- [ ] **Step 3: Verify all declared names resolve**

Run:
```bash
python -c "import npl.descriptors as d; [getattr(d,n) for n in d.__all__]; import npl.optimization as o; [getattr(o,n) for n in o.__all__]; print('exports ok')"
```
Expected: `exports ok` (raises `AttributeError` if any `__all__` entry is unbound).

### Task 2.3: Create the runnable flagship example (TOP energy evaluation)

**Files:**
- Create: `examples/top_energy_evaluation.py`

- [ ] **Step 1: Write the example script**

```python
"""Flagship example: evaluate the topological (TOP) energy of a bimetallic
nanoparticle using pretrained coefficients.

Builds a 201-atom truncated-octahedral Pt151Cu50 particle, computes its
extended topological feature vector, and evaluates the energy with a
TOPCalculator loaded from pretrained coefficients shipped with npl.
"""
from npl.core import Nanoparticle
from npl.calculators import TOPCalculator
from npl.descriptors import ExtendedTopologicalFeaturesClassifier


def main():
    calculator = TOPCalculator(
        "ETOP",
        stoichiometry="Pt151Cu50",
        feature_classifier=ExtendedTopologicalFeaturesClassifier,
    )
    feature_classifier = calculator.get_feature_classifier()

    particle = Nanoparticle()
    particle.truncated_octahedron(7, 2, {"Pt": 151, "Cu": 50})

    feature_classifier.compute_feature_vector(particle)
    energy = calculator.compute_energy(particle)

    print(f"TOP energy of Pt151Cu50: {float(energy):.4f}")
    return float(energy)


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run it end-to-end**

Run: `python examples/top_energy_evaluation.py`
Expected: prints `TOP energy of Pt151Cu50: <number>` and exits 0. If `compute_feature_vector` stores under a key other than `ETOP`, read `ExtendedTopologicalFeaturesClassifier.get_feature_key()` and pass that string as the first `TOPCalculator` arg so the keys match; re-run until it prints a finite number.

### Task 2.4: De-reference Monte Carlo in README and docs

**Files:**
- Modify: `README.md`
- Delete: `docs/npl.monte_carlo.rst`, `docs/npl.monte_carlo.ensembles.rst`
- Modify: `docs/npl.rst`, `docs/modules.rst` (remove `npl.monte_carlo` toctree entries if present)

- [ ] **Step 1: Replace the README "Monte Carlo Run Example" section**

In `README.md`, replace the entire "### Monte Carlo Run Example" block (the prose + the python code fence using `npl.monte_carlo`) with:

```markdown
### TOP Energy Evaluation Example

This example evaluates the topological (TOP) energy of a bimetallic
nanoparticle using pretrained coefficients. A 201-atom truncated-octahedral
Pt151Cu50 particle is built, its extended topological feature vector is
computed, and the energy is evaluated with a `TOPCalculator`.

​```python
from npl.core import Nanoparticle
from npl.calculators import TOPCalculator
from npl.descriptors import ExtendedTopologicalFeaturesClassifier

calculator = TOPCalculator('ETOP', stoichiometry='Pt151Cu50',
                           feature_classifier=ExtendedTopologicalFeaturesClassifier)
feature_classifier = calculator.get_feature_classifier()

particle = Nanoparticle()
particle.truncated_octahedron(7, 2, {'Pt': 151, 'Cu': 50})

feature_classifier.compute_feature_vector(particle)
energy = calculator.compute_energy(particle)
print(f"TOP energy: {float(energy):.4f}")
​```

> Looking for Monte Carlo / grand-canonical sampling? That functionality now
> lives in the companion package [mcpy](https://github.com/farrisric/mcpy).
```

(Remove the backtick zero-width markers when editing; they are only to show the fence here.) Also update the "About NPL" paragraph if it implies MC is part of npl.

- [ ] **Step 2: Remove MC doc pages and toctree references**

Run: `git rm docs/npl.monte_carlo.rst docs/npl.monte_carlo.ensembles.rst`
Then edit `docs/npl.rst` and `docs/modules.rst`: delete any line referencing `npl.monte_carlo` from their `toctree` lists.

- [ ] **Step 3: Verify no MC references remain in README/docs (text)**

Run: `grep -rn "monte_carlo\|run_monte_carlo\|mc_run" README.md docs/*.rst`
Expected: no output (notebook `.ipynb`/`docs/examples` handled next).

### Task 2.5: Migrate or retire the MC-based example notebooks

The notebooks `global_optimization.ipynb`, `multimet_go.ipynb`, `loading_pretrain_top.ipynb`, and `adsorbate_optimization.ipynb` import `npl.monte_carlo` (and one imports a non-existent `npl.monte_carlo.ensembles.canonical_ensemble`, so it is already broken). They are research notebooks; rewriting them to fully execute is out of scope. Retire them as npl examples and point users to mcpy, keeping only the runnable `top_energy_evaluation.py` plus any notebook that does not touch MC.

**Files:**
- Delete: the four MC-referencing notebooks under `examples/` and their `docs/examples/*.rst` mirrors
- Modify: `docs/examples.rst`, `docs/examples/` index entries

- [ ] **Step 1: Identify which example files reference MC**

Run: `grep -rln "monte_carlo\|run_monte_carlo\|mc_run\|canonical_ensemble" examples/ docs/examples/`
Record the list.

- [ ] **Step 2: Remove the MC-referencing examples and their rst mirrors**

For each file from Step 1, run `git rm <file>`. Also `git rm examples/npl.log` (a stray log artifact) and `git rm examples/GO_tab_data.ipynb` only if it references MC (check the grep output first; if it does not, leave it).

- [ ] **Step 3: Add a docs pointer to the new example and to mcpy**

In `docs/examples.rst`, replace the removed entries with a reference to `top_energy_evaluation` and a note: "Monte Carlo / GCMC examples are maintained in the companion package mcpy (https://github.com/farrisric/mcpy)."

- [ ] **Step 4: Verify the codebase is MC-free**

Run: `grep -rn "monte_carlo\|run_monte_carlo\|mc_run\|canonical_ensemble" npl/ examples/ docs/ README.md`
Expected: no output.

- [ ] **Step 5: Commit Phase 2**

```bash
git add -A
git commit -m "refactor: remove npl.monte_carlo, point users to mcpy

Delete the migrated-out Monte Carlo package and every example/doc that
referenced it, fix broken __all__ exports in descriptors/optimization,
and add a runnable TOP energy-evaluation flagship example."
```

---

# Phase 3 — Dead-code sweep, latent-bug fixes, operator clean-up

**Exit criterion:** `flake8` reports no errors; no unused public exports; the basin-hopping `flip` bug is fixed; both BH+TOP and BH+ACT paths are reachable via `run_basin_hopping(..., model=...)`.

### Task 3.1: Remove confirmed dead code

All verified unreferenced (`grep -rn <name> npl examples docs` shows only the definition):
- `compute_coefficients_for_shape_optimization` in `npl/calculators/energy_calculator.py` (never called)
- `LateralInteractionCalculator` in `energy_calculator.py` (only used by the now-deleted `npl/monte_carlo/monte_carlo.py`)
- `DipoleMomentCalculator` in `energy_calculator.py` (unused; a separate `DipoleMomentCalculator` exists in `global_feature_classifier.py` and is the one referenced elsewhere)
- The commented-out `brr_energy` lines in `BayesianRRCalculator.compute_energy` (`energy_calculator.py:330-332`)
- `data/params.json` (unused; `TOPCalculator` reads `npl/calculators/parameters.py:top_parameters`)

**DO NOT DELETE (user-confirmed, load-bearing despite no static importer):**
`npl/optimization/local_optimization/garbage_exchange_operator.py` is the
exchange operator for the **basin-hopping + ACT** workflow (the BH + TOP path
uses `guided_exchange_operator.py`). It shows no grep hits only because the
BH+ACT path is not wired into the runner. It is **renamed and wired** in
Tasks 3.3–3.4 and exercised in Phase 4 Task 4.5 — not removed.

**Files:**
- Delete: `data/params.json`
- Modify: `npl/calculators/energy_calculator.py`

- [ ] **Step 1: Re-verify each suspect is unreferenced**

Run:
```bash
for n in compute_coefficients_for_shape_optimization LateralInteractionCalculator; do
  echo "== $n =="; grep -rn "$n" npl/ examples/ docs/ | grep -v "energy_calculator.py:.*def ";
done
python -c "import json; json.load(open('data/params.json')); print('params.json is valid json but unused')"
```
Expected: each `==` block prints nothing (only definitions exist). If any prints a real use, keep that item and note it. (Do **not** add `garbage_exchange` to this loop — it is intentionally kept.)

- [ ] **Step 2: Delete the unused data file**

Run: `git rm data/params.json`

- [ ] **Step 3: Remove dead classes/function/comments from `energy_calculator.py`**

Delete the `compute_coefficients_for_shape_optimization` function, the `LateralInteractionCalculator` class, the `DipoleMomentCalculator` class, and the three commented-out `# brr_energy = ...` lines inside `BayesianRRCalculator.compute_energy`. Remove any imports left unused by these deletions (e.g. `copy` only if no longer used — verify with grep before removing).

- [ ] **Step 4: Verify calculators still import and the example still runs**

Run:
```bash
python -c "import npl.calculators; from npl.calculators import EMTCalculator, BayesianRRCalculator, TOPCalculator, EnergyCalculator; print('ok')"
python examples/top_energy_evaluation.py
```
Expected: `ok`, then a finite TOP energy.

### Task 3.2: Fix the basin-hopping `flip` expression bug

`npl/optimization/basin_hopping/basin_hopping.py:38-39` reads:
```python
flip = exchange_operator.symbol1_exchange_energies[index1]
+ exchange_operator.symbol2_exchange_energies[index2]
```
The second line is a separate no-op statement (unary `+`), so `flip` silently drops the second term.

**Files:**
- Modify: `npl/optimization/basin_hopping/basin_hopping.py`

- [ ] **Step 1: Join the expression into one statement**

Replace those two lines with:
```python
            flip = (exchange_operator.symbol1_exchange_energies[index1]
                    + exchange_operator.symbol2_exchange_energies[index2])
```

- [ ] **Step 2: Verify the module imports**

Run: `python -c "from npl.optimization.basin_hopping.basin_hopping import run_basin_hopping; print('ok')"`
Expected: `ok`.

### Task 3.3: Rename the ACT exchange operator and remove its dead commented code

`garbage_exchange_operator.py` is the basin-hopping + ACT exchange operator (user-confirmed, kept). Its class is confusingly named `GuidedExchangeOperator` — identical to the TOP operator's class in `guided_exchange_operator.py`. Rename the file and class for the JOSS release, and drop the dead commented-out second `update()` (lines ~150–204). Do **not** alter the live algorithm.

**Files:**
- Rename: `npl/optimization/local_optimization/garbage_exchange_operator.py` → `act_exchange_operator.py`
- Modify: the renamed file (class name + remove dead comment block)

- [ ] **Step 1: Rename the file**

Run: `git mv npl/optimization/local_optimization/garbage_exchange_operator.py npl/optimization/local_optimization/act_exchange_operator.py`

- [ ] **Step 2: Rename the class and delete the dead commented `update()`**

In `act_exchange_operator.py`, change `class GuidedExchangeOperator:` to `class ACTExchangeOperator:`. Delete the entire trailing commented-out block (the `# def update(self, particle, indices, exchange_indices):` section through the end of the file). Leave the live `__init__`, `env_from_feature`, `guided_exchange`, `basin_hop_step`, `bind_particle`, and the active `update` methods unchanged.

- [ ] **Step 3: Verify it imports under the new name**

Run: `python -c "from npl.optimization.local_optimization.act_exchange_operator import ACTExchangeOperator; print('ok')"`
Expected: `ok`. Also confirm the old name is gone: `grep -rn "garbage_exchange\|GuidedExchangeOperator" npl/optimization/local_optimization/act_exchange_operator.py` → no output.

### Task 3.4: Parameterize the runner to select TOP vs ACT operator

`setup_local_optimization` currently hardcodes `GuidedExchangeOperator` (TOP). Add a `model` selector so both basin-hopping paths are first-class.

**Files:**
- Modify: `npl/optimization/local_optimization/local_optimization.py`
- Modify: `npl/optimization/basin_hopping/basin_hopping.py`

- [ ] **Step 1: Add an operator registry and `model` parameter in `local_optimization.py`**

At the top of `npl/optimization/local_optimization/local_optimization.py`, add the ACT import and a registry next to the existing imports:

```python
from npl.optimization.local_optimization.guided_exchange_operator import GuidedExchangeOperator
from npl.optimization.local_optimization.act_exchange_operator import ACTExchangeOperator

EXCHANGE_OPERATORS = {"TOP": GuidedExchangeOperator, "ACT": ACTExchangeOperator}
```

Change `setup_local_optimization`'s signature to accept `model="TOP"` and select the operator class from the registry:

```python
def setup_local_optimization(start_particle, energy_calculator, environment_energies,
                             local_feature_classifier=None, model="TOP"):
    ...
    try:
        exchange_operator_cls = EXCHANGE_OPERATORS[model]
    except KeyError:
        raise ValueError(
            f"Unknown model {model!r}; choose from {sorted(EXCHANGE_OPERATORS)}")
    exchange_operator = exchange_operator_cls(environment_energies, feature_key)
    exchange_operator.bind_particle(start_particle)
    ...
```

(Keep the rest of the function body identical; only the operator instantiation line changes from the hardcoded class to `exchange_operator_cls`.)

- [ ] **Step 2: Thread `model` through `local_optimization` and `run_basin_hopping`**

In `local_optimization.py`, change `def local_optimization(start_particle, energy_calculator, environment_energies, local_feature_classifier=None):` to add `model="TOP"`, and pass `model=model` into its `setup_local_optimization(...)` call.

In `npl/optimization/basin_hopping/basin_hopping.py`, change `run_basin_hopping`'s signature to add `model="TOP"` (after `local_feature_classifier=None`) and pass `model=model` into its `setup_local_optimization(...)` call.

- [ ] **Step 3: Verify both operators are selectable**

Run:
```bash
python -c "from npl.optimization.local_optimization.local_optimization import EXCHANGE_OPERATORS; print(sorted(EXCHANGE_OPERATORS))"
python -c "import inspect; from npl.optimization.basin_hopping.basin_hopping import run_basin_hopping; print('model' in inspect.signature(run_basin_hopping).parameters)"
```
Expected: `['ACT', 'TOP']` then `True`.

### Task 3.5: Rename the misspelled `genetic_algoritm` package

**Files:**
- Rename: `npl/optimization/genetic_algoritm/` → `npl/optimization/genetic_algorithm/`
- Modify: internal imports in `SingleParticleGA.py`; `docs/npl.optimization.rst`, `docs/npl.optimization.genetic_algoritm.rst`

- [ ] **Step 1: Rename the directory**

Run: `git mv npl/optimization/genetic_algoritm npl/optimization/genetic_algorithm`

- [ ] **Step 2: Update internal references**

Run: `grep -rln "genetic_algoritm" npl/ docs/` and in each hit replace `genetic_algoritm` with `genetic_algorithm` (notably `npl/optimization/genetic_algorithm/SingleParticleGA.py` lines 2–3, and `docs/npl.optimization.rst`). Rename `docs/npl.optimization.genetic_algoritm.rst` → `docs/npl.optimization.genetic_algorithm.rst` with `git mv` and fix the `automodule` paths inside it.

- [ ] **Step 3: Verify**

Run: `grep -rn "genetic_algoritm" npl/ docs/` → expect no output. Then `python -c "import npl.optimization.genetic_algorithm.SingleParticleGA"` (if it has importable deps) or at minimum `python -c "import npl.optimization; print('ok')"`.
Expected: `ok`, no stale spelling.

### Task 3.6: flake8 clean + commit

- [ ] **Step 1: Run the linter**

Run: `flake8 npl/`
Expected: no output (config in `.flake8`: max line 100). Fix any unused-import / line-length issues introduced by the deletions. Note `act_exchange_operator.py` and `guided_exchange_operator.py` have pre-existing long lines; if flake8 now flags them, wrap to ≤100 cols without changing logic.

- [ ] **Step 2: Commit**

```bash
git add -A
git commit -m "refactor: remove dead code, fix latent bugs, clean up exchange operators

Drop unused calculators (LateralInteraction, DipoleMoment in
energy_calculator), compute_coefficients_for_shape_optimization, and
unused data/params.json; fix the basin-hopping flip expression bug;
rename genetic_algoritm -> genetic_algorithm; rename the ACT exchange
operator (garbage_exchange_operator -> act_exchange_operator,
GuidedExchangeOperator -> ACTExchangeOperator) and parameterize
run_basin_hopping/local_optimization with a TOP/ACT model selector."
```

---

# Phase 4 — Test suite + CI

**Exit criterion:** `pytest` passes locally; the CI workflow runs pytest green on a Python matrix.

### Task 4.1: Create the test directory and import test

**Files:**
- Delete: `test/test_core.py` (empty placeholder), and the `test/` dir
- Create: `tests/conftest.py`, `tests/test_imports.py`

- [ ] **Step 1: Remove the empty legacy test dir**

Run: `git rm -r test`

- [ ] **Step 2: Write `tests/conftest.py` with shared fixtures**

```python
import pytest
from npl.core import Nanoparticle


@pytest.fixture
def ptcu_particle():
    """A 201-atom truncated-octahedral Pt151Cu50 nanoparticle."""
    particle = Nanoparticle()
    particle.truncated_octahedron(7, 2, {"Pt": 151, "Cu": 50})
    return particle
```

- [ ] **Step 3: Write `tests/test_imports.py`**

```python
import importlib
import pytest

SUBMODULES = [
    "npl",
    "npl.core",
    "npl.descriptors",
    "npl.calculators",
    "npl.optimization",
    "npl.visualize",
    "npl.utils",
]


@pytest.mark.parametrize("module_name", SUBMODULES)
def test_submodule_imports(module_name):
    importlib.import_module(module_name)


def test_no_monte_carlo_module():
    with pytest.raises(ModuleNotFoundError):
        importlib.import_module("npl.monte_carlo")


def test_declared_exports_resolve():
    import npl.descriptors as d
    import npl.optimization as o
    for name in d.__all__:
        assert getattr(d, name) is not None
    for name in o.__all__:
        assert getattr(o, name) is not None
```

- [ ] **Step 4: Run**

Run: `pytest tests/test_imports.py -v`
Expected: all PASS.

### Task 4.2: Core tests

**Files:**
- Create: `tests/test_core.py`

- [ ] **Step 1: Write the tests**

```python
def test_truncated_octahedron_atom_count(ptcu_particle):
    assert ptcu_particle.atoms.get_n_atoms() == 201


def test_truncated_octahedron_stoichiometry(ptcu_particle):
    stoich = ptcu_particle.get_stoichiometry()
    assert stoich["Pt"] == 151
    assert stoich["Cu"] == 50


def test_neighbor_list_is_built(ptcu_particle):
    indices = list(ptcu_particle.get_indices())
    assert len(indices) == 201
    assert all(len(ptcu_particle.neighbor_list[i]) > 0 for i in indices)


def test_geometrical_data_has_positions_and_neighbors(ptcu_particle):
    data = ptcu_particle.get_geometrical_data()
    assert "positions" in data and "neighbor_list" in data
    assert len(data["positions"]) == 201
```

- [ ] **Step 2: Run and adjust**

Run: `pytest tests/test_core.py -v`
Expected: PASS. If a method name differs from the assumed public API (`get_n_atoms`, `get_stoichiometry`, `get_indices`, `get_geometrical_data`), confirm the real name with `python -c "import npl.core.base_nanoparticle as b; print([m for m in dir(b.BaseNanoparticle) if not m.startswith('_')])"` and update the test to the actual method.

### Task 4.3: Descriptor tests

**Files:**
- Create: `tests/test_descriptors.py`

- [ ] **Step 1: Write the tests**

```python
import numpy as np
from npl.descriptors import (
    ExtendedTopologicalFeaturesClassifier,
    NeighborCountingEnvironmentCalculator,
)


def test_neighbor_counting_environment(ptcu_particle):
    calc = NeighborCountingEnvironmentCalculator(["Pt", "Cu"])
    env = calc.predict_local_environment(ptcu_particle, 0)
    assert len(env) == 2
    # the two counts sum to the coordination number of atom 0
    assert env.sum() == len(ptcu_particle.neighbor_list[0])


def test_extended_topological_feature_vector(ptcu_particle):
    classifier = ExtendedTopologicalFeaturesClassifier(["Pt", "Cu"])
    classifier.compute_feature_vector(ptcu_particle)
    key = classifier.get_feature_key()
    fv = ptcu_particle.get_feature_vector(key)
    assert fv is not None
    assert np.asarray(fv).ndim == 1
    assert len(fv) == len(classifier.get_feature_labels())
```

- [ ] **Step 2: Run and adjust**

Run: `pytest tests/test_descriptors.py -v`
Expected: PASS. If `ExtendedTopologicalFeaturesClassifier.__init__` needs a different symbols argument shape, inspect with `python -c "import inspect, npl.descriptors as d; print(inspect.signature(d.ExtendedTopologicalFeaturesClassifier.__init__))"` and adjust.

### Task 4.4: Calculator tests

**Files:**
- Create: `tests/test_calculators.py`

- [ ] **Step 1: Write the tests**

```python
import numpy as np
from npl.core import Nanoparticle
from npl.calculators import TOPCalculator, EMTCalculator
from npl.descriptors import ExtendedTopologicalFeaturesClassifier


def test_top_calculator_loads_and_evaluates(ptcu_particle):
    calc = TOPCalculator(
        "ETOP",
        stoichiometry="Pt151Cu50",
        feature_classifier=ExtendedTopologicalFeaturesClassifier,
    )
    fc = calc.get_feature_classifier()
    fc.compute_feature_vector(ptcu_particle)
    energy = calc.compute_energy(ptcu_particle)
    assert np.isfinite(float(energy))


def test_emt_calculator_sets_energy():
    particle = Nanoparticle()
    particle.truncated_octahedron(7, 2, {"Au": 101, "Ag": 100})
    EMTCalculator().compute_energy(particle)
    assert np.isfinite(particle.get_energy("EMT"))
```

- [ ] **Step 2: Run and adjust**

Run: `pytest tests/test_calculators.py -v`
Expected: PASS. If the TOP feature key differs from `"ETOP"`, set the first arg of `TOPCalculator` to `fc.get_feature_key()` after constructing `fc` (mirror the fix from Task 2.3 Step 2). EMT supports Au/Ag, so the EMT test needs no relaxation (default `relax_atoms=False`).

### Task 4.5: Optimization integration test

`local_optimization(start_particle, energy_calculator, environment_energies, local_feature_classifier=None)` runs a guided-exchange local optimization and returns `[particle, accepted_energies]`. `environment_energies` is the per-environment coefficient array from `compute_coefficients_for_linear_topological_model`.

**Files:**
- Create: `tests/test_optimization.py`

- [ ] **Step 1: Write a smoke test that the local optimizer runs and does not increase energy**

```python
import numpy as np
from npl.core import Nanoparticle
from npl.calculators import (
    TOPCalculator,
    compute_coefficients_for_linear_topological_model,
)
from npl.descriptors import ExtendedTopologicalFeaturesClassifier
from npl.optimization import local_optimization


def test_local_optimization_runs_and_does_not_worsen():
    symbols = ["Pt", "Cu"]
    calc = TOPCalculator(
        "ETOP",
        stoichiometry="Pt151Cu50",
        feature_classifier=ExtendedTopologicalFeaturesClassifier,
    )
    fc = calc.get_feature_classifier()

    particle = Nanoparticle()
    particle.truncated_octahedron(7, 2, {"Pt": 151, "Cu": 50})
    fc.compute_feature_vector(particle)

    coefficients = calc.get_coefficients()
    environment_energies, _ = compute_coefficients_for_linear_topological_model(
        coefficients, symbols, n_atoms=201
    )

    result_particle, accepted = local_optimization(particle, calc, environment_energies)

    assert len(accepted) >= 1
    # accepted energies are (energy, step) tuples; the last is <= the first
    assert accepted[-1][0] <= accepted[0][0] + 1e-6
```

- [ ] **Step 2: Run and adjust**

Run: `pytest tests/test_optimization.py -v`
Expected: PASS. This test couples several subsystems; if the `environment_energies` shape expected by `GuidedExchangeOperator` differs from the output of `compute_coefficients_for_linear_topological_model`, read `guided_exchange_operator.py` and `setup_local_optimization` to determine the exact array the operator indexes, and construct `environment_energies` accordingly. If the coupling proves impractical to satisfy, downgrade this to a "runs without raising" assertion (drop the monotonicity check) but keep the test executing the real `local_optimization` path.

- [ ] **Step 3: Add a test for the basin-hopping + ACT path**

After Phase 3 Task 3.4, `run_basin_hopping(..., model="ACT")` selects the `ACTExchangeOperator`. Add a test that runs both basin-hopping models on the same particle with tiny iteration counts and returns a best particle:

```python
import pytest
from npl.optimization.basin_hopping.basin_hopping import run_basin_hopping


@pytest.mark.parametrize("model", ["TOP", "ACT"])
def test_basin_hopping_runs_for_both_models(model):
    symbols = ["Pt", "Cu"]
    calc = TOPCalculator(
        "ETOP",
        stoichiometry="Pt151Cu50",
        feature_classifier=ExtendedTopologicalFeaturesClassifier,
    )
    fc = calc.get_feature_classifier()
    particle = Nanoparticle()
    particle.truncated_octahedron(7, 2, {"Pt": 151, "Cu": 50})
    fc.compute_feature_vector(particle)

    environment_energies, _ = compute_coefficients_for_linear_topological_model(
        calc.get_coefficients(), symbols, n_atoms=201
    )

    best, lowest_energies, _ = run_basin_hopping(
        particle, calc, environment_energies,
        n_hopping_attempts=2, n_hops=2, model=model,
    )
    assert best is not None
    assert len(lowest_energies) >= 1
```

Expected: PASS for both `model` values. The ACT operator (`ACTExchangeOperator`) additionally reads per-atom integer features via `particle.get_coordination_atoms(...)`; if the default `TopologicalEnvironmentClassifier` features it consumes are not the integer coordination-type features ACT expects, the ACT case will raise. If so, set the correct local feature classifier for ACT (see open item #4) by passing `local_feature_classifier=` and re-run; if the ACT inputs cannot be resolved without further domain input, mark the ACT parametrization `xfail` with a reason referencing open item #4 rather than deleting it.

### Task 4.6: Enable pytest in CI + commit

**Files:**
- Modify: `.github/workflows/test.yml`

- [ ] **Step 1: Rewrite the CI workflow to run pytest on a matrix**

Replace the contents of `.github/workflows/test.yml` with:

```yaml
name: CI

on:
  push:
    branches: [main]
  pull_request:
    branches: [main]

jobs:
  test:
    runs-on: ubuntu-latest
    strategy:
      fail-fast: false
      matrix:
        python-version: ["3.10", "3.11", "3.12", "3.13"]

    steps:
      - uses: actions/checkout@v4
      - name: Set up Python ${{ matrix.python-version }}
        uses: actions/setup-python@v5
        with:
          python-version: ${{ matrix.python-version }}
      - name: Install package with test extras
        run: |
          python -m pip install --upgrade pip
          pip install -e ".[test]"
      - name: Run tests
        run: pytest -v
```

(The `[test]` extra is defined in Phase 5's `pyproject.toml`. Until then, run tests locally.)

- [ ] **Step 2: Run the whole suite locally**

Run: `pytest -v`
Expected: all tests PASS.

- [ ] **Step 3: Commit**

```bash
git add -A
git commit -m "test: add pytest suite for core, descriptors, calculators, optimization

Replace the empty placeholder test with real coverage and enable pytest
in CI across Python 3.10-3.13."
```

---

# Phase 5 — Packaging + metadata

**Exit criterion:** `python -m build` produces an sdist and wheel; installing the wheel into a clean venv imports `npl` and reports version `1.1.0`; `pip install -e ".[test]"` works.

### Task 5.1: Single-source the version

**Files:**
- Modify: `npl/__init__.py`

- [ ] **Step 1: Set the canonical version to 1.1.0**

In `npl/__init__.py`, change `__version__ = "1.0.0"` to `__version__ = "1.1.0"`. (1.1.0 reflects the MC removal + modernization; the MC removal is a breaking change to the public surface, so if the user prefers, this can be `2.0.0` — confirm before tagging a release.)

### Task 5.2: Replace `setup.py` with `pyproject.toml`

**Files:**
- Create: `pyproject.toml`
- Delete: `setup.py`

- [ ] **Step 1: Write `pyproject.toml`**

```toml
[project]
name = "npl"
description = "NanoParticleLibrary: structural and chemical-ordering analysis and optimization of nanoparticles"
readme = "README.md"
requires-python = ">=3.10"
license = { text = "MIT" }
authors = [{ name = "Riccardo Farris", email = "rfarris@ub.edu" }]
dynamic = ["version"]
dependencies = [
    "ase>=3.23",
    "numpy>=1.26",
    "scipy>=1.11",
    "scikit-learn>=1.5",
    "sortedcontainers>=2.4",
]

[project.optional-dependencies]
test = ["pytest>=7.0"]

[project.urls]
Homepage = "https://github.com/farrisric/npl"
Documentation = "https://nplib.readthedocs.io"
Issues = "https://github.com/farrisric/npl/issues"

[build-system]
requires = ["setuptools>=64", "wheel"]
build-backend = "setuptools.build_meta"

[tool.setuptools.packages.find]
include = ["npl*"]

[tool.setuptools.dynamic]
version = { attr = "npl.__version__" }

[tool.pytest.ini_options]
testpaths = ["tests"]
```

Note: the `acat` dependency was in `setup.py`/`requirements.txt` but is not imported anywhere in `npl/` (verify: `grep -rn "import acat\|from acat" npl/`). Omit it unless that grep finds a real import; if it does, add `"acat>=1.7.1"` back to `dependencies`.

- [ ] **Step 2: Remove setup.py**

Run: `git rm setup.py`

- [ ] **Step 3: Reinstall and verify version is single-sourced**

Run:
```bash
pip install -e ".[test]"
python -c "import npl; print(npl.__version__)"
pip show npl | grep -i version
```
Expected: both report `1.1.0`.

### Task 5.3: Update the publish workflow

**Files:**
- Modify: `.github/workflows/publish.yml`

- [ ] **Step 1: Switch the build step to `python -m build`**

In `.github/workflows/publish.yml`, replace the "Install dependencies" + "Build package" steps with:

```yaml
      - name: Install build tooling
        run: |
          python -m pip install --upgrade pip
          pip install build
      - name: Build package
        run: python -m build
```

Bump `python-version` in that workflow to `"3.11"`.

- [ ] **Step 2: Verify a clean build**

Run:
```bash
pip install build
python -m build
ls dist/
```
Expected: `dist/` contains `npl-1.1.0.tar.gz` and `npl-1.1.0-py3-none-any.whl`.

- [ ] **Step 3: Update README/requirements for the new floor**

In `README.md`, change "Python 3.9+" to "Python 3.10+". Update `requirements.txt` to match the `pyproject.toml` floors (or delete it and rely on `pip install -e ".[test]"`; if deleting, update the README/CI references to it — note CI already uses `pip install -e ".[test]"`).

- [ ] **Step 4: Commit**

```bash
git add -A
git commit -m "build: migrate to pyproject.toml and single-source version 1.1.0

Replace setup.py with PEP 621 pyproject.toml, align dependency floors
with the supported stack, drop the unused acat dependency, and build the
publish workflow with python -m build."
```

---

# Phase 6 — JOSS paper + repo essentials

**Exit criterion:** `paper/paper.md` has valid JOSS YAML front-matter, every `paper.bib` key it cites resolves, and the repo has `CONTRIBUTING.md` plus an accurate `CITATION.cff`.

### Task 6.1: Write `paper/paper.bib`

**Files:**
- Create: `paper/paper.bib`

- [ ] **Step 1: Write the bibliography**

```bibtex
@article{farris2024ordering,
  author  = {Farris, Riccardo and Neyman, Konstantin M. and Bruix, Albert},
  title   = {Determining the chemical ordering in nanoalloys by considering atomic coordination types},
  journal = {The Journal of Chemical Physics},
  volume  = {161},
  number  = {13},
  pages   = {134114},
  year    = {2024},
  doi     = {10.1063/5.0214377}
}

@article{farris2024zr,
  author  = {Farris, Riccardo and Merinov, Boris V. and Bruix, Albert and Neyman, Konstantin M.},
  title   = {Effects of {Zr} dopants on properties of {PtNi} nanoparticles for {ORR} catalysis: A {DFT} modeling},
  journal = {The Journal of Chemical Physics},
  volume  = {160},
  number  = {12},
  pages   = {124706},
  year    = {2024},
  doi     = {10.1063/5.0193848}
}

@article{neumann2021,
  author  = {Neumann, Felix and Margraf, Johannes T. and Reuter, Karsten and Bruix, Albert},
  title   = {Interplay between shape and composition in bimetallic nanoparticles revealed by an efficient optimal-exchange optimization algorithm},
  journal = {ChemRxiv},
  year    = {2021},
  doi     = {10.26434/chemrxiv-2021-26ztp}
}

@article{ase2017,
  author  = {Larsen, Ask Hjorth and Mortensen, Jens J{\o}rgen and Blomqvist, Jakob and others},
  title   = {The Atomic Simulation Environment---a {Python} Library for Working with Atoms},
  journal = {Journal of Physics: Condensed Matter},
  volume  = {29},
  number  = {27},
  pages   = {273002},
  year    = {2017},
  doi     = {10.1088/1361-648X/aa680e}
}

@article{scikit-learn,
  author  = {Pedregosa, F. and Varoquaux, G. and Gramfort, A. and others},
  title   = {Scikit-learn: Machine Learning in {Python}},
  journal = {Journal of Machine Learning Research},
  volume  = {12},
  pages   = {2825--2830},
  year    = {2011}
}

@misc{mcpy,
  author       = {Farris, Riccardo and Telari, Emanuele and Bruix, Albert},
  title        = {{mcpy}: Grand Canonical Monte Carlo for atomistic systems with machine-learning potentials},
  year         = {2026},
  howpublished = {\url{https://github.com/farrisric/mcpy}}
}
```

### Task 6.2: Write `paper/paper.md`

**Files:**
- Create: `paper/paper.md`

- [ ] **Step 1: Write the paper**

```markdown
---
title: 'NanoParticleLibrary: topological descriptors and chemical-ordering optimization for bimetallic nanoparticles'
tags:
  - Python
  - computational chemistry
  - materials science
  - nanoparticles
  - nanoalloys
  - chemical ordering
  - heterogeneous catalysis
authors:
  - name: Riccardo Farris
    orcid: 0000-0002-8157-6786
    affiliation: "1, 2"
  - name: Albert Bruix
    orcid: 0000-0003-2585-5542
    affiliation: 1
affiliations:
  - name: Institut de Química Teòrica i Computacional (IQTCUB), Universitat de Barcelona, Barcelona, Spain
    index: 1
  - name: LEITAT Technological Center, Terrassa, Spain
    index: 2
date: 29 May 2026
bibliography: paper.bib
---

# Summary

`NanoParticleLibrary` (`npl`) is a Python package for building bimetallic
nanoparticles, describing them through coordination-based topological
descriptors, and optimizing their chemical ordering. It is built on the Atomic
Simulation Environment (ASE) [@ase2017], representing each particle through a
neighbor list and per-atom local environments from which feature vectors are
computed. These feature vectors feed a linear topological energy model whose
coefficients are fitted with scikit-learn [@scikit-learn] regressors (Bayesian
ridge regression) against reference energies, or loaded from pretrained
parameters shipped with the package. Given an energy model, `npl` searches for
low-energy chemical orderings with guided local optimization, basin hopping,
and genetic-algorithm strategies. The topological descriptors and the
optimal-exchange ideas behind these methods were introduced in
[@neumann2021; @farris2024ordering] and applied to catalytic nanoparticles in
[@farris2024zr].

# Statement of need

The properties of bimetallic nanoparticles depend not only on their size and
shape but strongly on their chemical ordering — how the two elements are
arranged among the available sites. Exhaustively enumerating orderings is
combinatorially impossible, and evaluating each candidate with density
functional theory is prohibitive. `npl` addresses this by combining cheap,
physically motivated topological descriptors with fitted energy models and
dedicated ordering optimizers, so that researchers can predict and explore the
most stable arrangements of nanoalloys at low cost. It is aimed at
computational chemists and materials scientists studying nanoalloys for
heterogeneous catalysis, where chemical ordering controls activity and
stability.

`npl` focuses on structure and chemical-ordering optimization with a fixed
number of atoms. Sampling at controlled chemical potential — grand canonical
Monte Carlo with machine-learning potentials — is provided by the companion
package `mcpy` [@mcpy]; the two packages share an ASE-native design and are
intended to be used together.

# Software design

`npl` is organized around four extensible abstractions connected by string
keys on a central particle object. A `Nanoparticle` holds atoms, a neighbor
list, and dictionaries — keyed by name — of energies, feature vectors, and
local environments, so that several descriptors and energy models can coexist
on one particle and geometry can be reused across different orderings.
*Descriptors* compute global feature vectors (topological and
coordination-based classifiers) or per-atom local environments. *Calculators*
implement a common `compute_energy(particle)` interface and write into the
particle's energy dictionary; they include an ASE EMT calculator, scikit-learn
regressors over feature vectors, and the linear topological model with
pretrained coefficients. *Optimizers* search chemical-ordering space using
guided exchanges, basin hopping, and genetic algorithms. This separation lets
users swap energy models and search strategies without touching the particle
representation.

# Acknowledgements

The authors received no specific funding for this work.

# AI usage disclosure

Claude (Opus) assisted with code modernization, dead-code removal, test
scaffolding, packaging configuration, and editing of this manuscript and the
documentation. All AI-assisted outputs were reviewed, edited, and validated by
the human authors, who designed the code architecture and take full
responsibility for the accuracy and originality of the software and all
submitted materials.

# References
```

- [ ] **Step 2: Verify front-matter and citation keys**

Run:
```bash
python -c "import yaml,io,sys; t=open('paper/paper.md').read(); fm=t.split('---')[1]; d=yaml.safe_load(fm); assert d['title'] and d['authors'] and d['bibliography']=='paper.bib'; print('front-matter ok')"
grep -oE '@[a-zA-Z0-9_]+' paper/paper.md | tr -d '@' | sort -u
grep -oE '^@[a-z]+\{[a-zA-Z0-9_]+' paper/paper.bib | sed -E 's/^@[a-z]+\{//' | sort -u
```
Expected: `front-matter ok`, and every citation key from the first list appears in the second (bib) list. (`pyyaml` is available in the environment.)

### Task 6.3: Add CONTRIBUTING.md and verify CITATION.cff

**Files:**
- Create: `CONTRIBUTING.md`
- Modify: `CITATION.cff` (only if inaccurate)

- [ ] **Step 1: Write `CONTRIBUTING.md`** (mirrors mcpy's structure)

```markdown
# Contributing to npl

Thanks for your interest in improving `npl`. Contributions of all kinds are
welcome: bug reports, feature requests, documentation, and code.

## Reporting bugs and requesting features

Please open an issue on the
[GitHub issue tracker](https://github.com/farrisric/npl/issues). For bugs,
include what you expected, what happened, a minimal reproducing example, and
your OS, Python version, and `npl` version.

## Development setup

```sh
git clone https://github.com/farrisric/npl.git
cd npl
pip install -e ".[test]"
```

## Running the tests

```sh
pytest
```

CI runs the test suite on Python 3.10–3.13.

## Linting

```sh
flake8 npl/
```

## Pull requests

Keep pull requests focused, add tests for new behavior, and make sure
`pytest` and `flake8` pass before requesting review.
```

- [ ] **Step 2: Check `CITATION.cff` accuracy**

Run: `cat CITATION.cff` and confirm the version field reads `1.1.0` and author/title match the project. If `cffconvert` is available, run `cffconvert --validate`; otherwise validate as YAML: `python -c "import yaml; yaml.safe_load(open('CITATION.cff')); print('cff yaml ok')"`. Update the version field to `1.1.0` if it differs.

- [ ] **Step 3: Commit**

```bash
git add -A
git commit -m "docs: add JOSS paper draft, CONTRIBUTING, and citation metadata

Draft paper/paper.md and paper/paper.bib (cross-citing mcpy), add
CONTRIBUTING.md, and align CITATION.cff with version 1.1.0."
```

---

## Final verification (run after all phases)

- [ ] `pytest -v` → all green
- [ ] `flake8 npl/` → no output
- [ ] `python -m build` → sdist + wheel in `dist/`
- [ ] `grep -rn "monte_carlo\|np.int\|set_calculator\|sph_harm\|genetic_algoritm\|garbage_exchange" npl/ examples/ README.md` → no output
- [ ] `python examples/top_energy_evaluation.py` → finite energy
- [ ] paper front-matter + bib keys validated (Task 6.2 Step 2)

## Open items needing the user's confirmation during execution

1. **Paper author list / affiliations** — drafted as Farris + Bruix mirroring mcpy and npl's own reference papers; confirm the final author list, ORCIDs, and affiliations.
2. **AI usage disclosure** — included in `paper.md` per mcpy's precedent (disclosure, not authorship). Remove if the target venue/user prefers not to disclose.
3. **Version number** — `1.1.0` proposed; `2.0.0` is defensible since removing `npl.monte_carlo` is a breaking API change.
4. **ACT feature-classifier input** *(resolved at the API level; one physics detail remains)* — the runner is now parameterized (`model="TOP"|"ACT"`, Task 3.4) and the operator renamed to `ACTExchangeOperator` (Task 3.3). The remaining unknown is which local feature classifier produces the integer coordination-type atom features the `ACTExchangeOperator` consumes via `get_coordination_atoms`. If the default `TopologicalEnvironmentClassifier` is not it, the BH+ACT test (Task 4.5 Step 3) needs the correct classifier passed as `local_feature_classifier=`. Confirm the intended ACT feature classifier (likely a coordination-types classifier from `npl.descriptors`) during execution.
