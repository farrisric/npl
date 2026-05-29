# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

NanoParticleLibrary (`npl`) — a Python library for simulating and optimizing the structure and
chemical ordering of (primarily bimetallic) nanoparticles. It is built on top of ASE (Atomic
Simulation Environment) and targets Python 3.9+.

## Commands

```bash
pip install -e .                  # editable install
pip install -r requirements.txt   # pinned dev/runtime deps (acat, scikit-learn, nbconvert, ...)

flake8                            # lint; config in .flake8 (max line 100, single quotes, E203 ignored)

python setup.py sdist bdist_wheel # build distributions (PyPI publish runs on pushing a v* tag)

cd docs && make html              # build Sphinx docs into docs/_build/
```

**Tests**: there is effectively no test suite yet. `test/test_core.py` is empty and the `pytest`
step in `.github/workflows/test.yml` is commented out — CI only installs deps. Don't assume a
passing test run validates a change; verify behavior through the example notebooks instead.

**Examples are the real documentation.** The tutorials live as Jupyter notebooks in `examples/`
(`global_optimization.ipynb`, `train_top.ipynb`, `loading_pretrain_top.ipynb`, `multimet_go.ipynb`,
`adsorbate_optimization.ipynb`) and are mirrored as `.rst` under `docs/examples/`.

## Architecture

The standard workflow is a pipeline tied together by **string keys**:

```
build Nanoparticle  →  feature_classifier.compute_feature_vector(p)  →
energy_calculator.compute_energy(p)  →  optimization / monte_carlo
```

A `feature_key` identifies which feature vector a calculator reads; an `energy_key` identifies where
a calculator writes its result. Calculators assume the matching feature vector has already been
computed on the particle.

- **`npl.core`** — the data model. `BaseNanoparticle` is a *data holder*: it wraps an `AtomWrapper`
  (atoms/symbols/positions), a `NeighborList`, an `AdsorptionSiteList`, and dictionaries keyed by
  string for `energies`, `feature_vectors`, `local_environments`, and `atom_features`. The
  string-keyed dicts are the core design choice — multiple energy models and descriptors coexist on
  one particle. Geometry (positions + neighbor list) is deliberately separable from ordering
  (symbols), so one geometry can back many orderings (`get_geometrical_data`, `save_npl_format`).
  `Nanoparticle(BaseNanoparticle)` adds project-specific construction such as
  `truncated_octahedron(...)` plus adsorption-site handling. Keep generic logic in `BaseNanoparticle`
  and project-specific logic in `Nanoparticle`.

- **`npl.descriptors`** — feature computation, two families: *global* classifiers (e.g.
  `TopologicalFeatureClassifier`, `ExtendedTopologicalFeaturesClassifier`,
  `CoordinationFeatureClassifier`) produce a particle-level `feature_vector`; *local* classifiers
  and calculators (e.g. `NeighborCountingEnvironmentCalculator`, `TopologicalEnvironmentClassifier`)
  produce per-atom local environments. The "topological"/coordination-type descriptors are the basis
  of the linear energy model used throughout (see the cited papers in README).

- **`npl.calculators`** — energy models. All subclass `EnergyCalculator` and implement
  `compute_energy(particle)`, writing into `particle.energies[energy_key]`. Implementations:
  `EMTCalculator` (ASE EMT + optional BFGS relaxation), `BayesianRRCalculator` and `GPRCalculator`
  (ML regressors over feature vectors), `MixingEnergyCalculator`, and `TOPCalculator` (the linear
  topological model). `TOPCalculator` loads pre-trained coefficients **by stoichiometry string**
  from the `top_parameters` dict in `npl/calculators/parameters.py`, or from a pickled model via
  `model_paths`. `compute_coefficients_for_linear_topological_model` converts fitted linear/BRR
  coefficients into the per-coordination-type topological coefficients.

- **`npl.optimization`** — global-optimization drivers. `GOSearch` is the base, specialized by
  `MCSearch`, `GASearch`, `GuidedSearch` (each fits an energy expression then runs a search and can
  repeat runs via `run_multiple_simulations`). Subpackages: `basin_hopping`, `local_optimization`,
  and `genetic_algoritm` (note the spelling — that is the actual directory name).

- **`npl.monte_carlo`** — Metropolis MC over chemical ordering. `run_monte_carlo(temperature,
  max_steps, start_particle, energy_calculator, feature_classifier)` is the main entry point. It uses
  random exchange operators and **incremental feature updates** — after a swap only the neighborhood
  of the exchanged atoms is recomputed (`features_to_update`), so don't introduce full
  recomputation in the MC loop.

- **`npl.visualize`**, **`npl.utils`** — plotting helpers and math/geometry utilities.

## Gotchas

- **The atomistic Monte Carlo work moved to a separate repo, `mcpy`** (commit `6616c23`,
  https://github.com/farrisric/mcpy). The additional working directory
  `/home/energystorage/projects/mcpy/mcpy/calculators` is its home and contains newer calculators
  (MACE, ALIGNN/alchemi). New MC / ML-potential work generally belongs in `mcpy`, not here.
- `import npl` runs `setup_logging()`, which writes an `npl.log` file in the current working
  directory and logs to stdout at INFO. Several modules also call `logging.basicConfig` again.
- `data/params.json` exists but `TOPCalculator` reads its coefficients from the Python dict in
  `npl/calculators/parameters.py`, not from that JSON file.
- Version is inconsistent: `setup.py` declares `1.0.4` while `npl/__init__.py` sets `__version__ =
  "1.0.0"`.
- Some code uses APIs removed in modern NumPy (e.g. `np.int` in `npl/core/nanoparticle.py`); be aware
  of NumPy/scikit-learn version sensitivity when running older paths.
Behavioral guidelines to reduce common LLM coding mistakes. Merge with project-specific instructions as needed.

**Tradeoff:** These guidelines bias toward caution over speed. For trivial tasks, use judgment.

## 1. Think Before Coding

**Don't assume. Don't hide confusion. Surface tradeoffs.**

Before implementing:
- State your assumptions explicitly. If uncertain, ask.
- If multiple interpretations exist, present them - don't pick silently.
- If a simpler approach exists, say so. Push back when warranted.
- If something is unclear, stop. Name what's confusing. Ask.

## 2. Simplicity First

**Minimum code that solves the problem. Nothing speculative.**

- No features beyond what was asked.
- No abstractions for single-use code.
- No "flexibility" or "configurability" that wasn't requested.
- No error handling for impossible scenarios.
- If you write 200 lines and it could be 50, rewrite it.

Ask yourself: "Would a senior engineer say this is overcomplicated?" If yes, simplify.

## 3. Surgical Changes

**Touch only what you must. Clean up only your own mess.**

When editing existing code:
- Don't "improve" adjacent code, comments, or formatting.
- Don't refactor things that aren't broken.
- Match existing style, even if you'd do it differently.
- If you notice unrelated dead code, mention it - don't delete it.

When your changes create orphans:
- Remove imports/variables/functions that YOUR changes made unused.
- Don't remove pre-existing dead code unless asked.

The test: Every changed line should trace directly to the user's request.

## 4. Goal-Driven Execution

**Define success criteria. Loop until verified.**

Transform tasks into verifiable goals:
- "Add validation" → "Write tests for invalid inputs, then make them pass"
- "Fix the bug" → "Write a test that reproduces it, then make it pass"
- "Refactor X" → "Ensure tests pass before and after"

For multi-step tasks, state a brief plan:
```
1. [Step] → verify: [check]
2. [Step] → verify: [check]
3. [Step] → verify: [check]
```

Strong success criteria let you loop independently. Weak criteria ("make it work") require constant clarification.

---

**These guidelines are working if:** fewer unnecessary changes in diffs, fewer rewrites due to overcomplication, and clarifying questions come before implementation rather than after mistakes.