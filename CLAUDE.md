# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

NanoParticleLibrary (`npl`) — a Python library for simulating and optimizing the structure and
chemical ordering of (primarily bimetallic) nanoparticles. It is built on top of ASE (Atomic
Simulation Environment) and targets Python 3.10+.

## Commands

```bash
pip install -e ".[test]"          # editable install + test deps (pytest); deps live in pyproject.toml

pytest                            # run the test suite (config in pyproject.toml: testpaths=tests)
flake8                            # lint; config in .flake8 (max line 100, single quotes, E203 ignored)

python -m build                   # build sdist+wheel (PyPI publish runs on pushing a v* tag)

cd docs && make html              # build Sphinx docs into docs/_build/
```

**Tests**: a real `pytest` suite lives in `tests/` (imports, core, descriptors, calculators,
optimization) and CI (`.github/workflows/test.yml`) runs it on Python 3.10–3.13. Add/adjust tests
when changing behavior.

**Examples**: `examples/top_energy_evaluation.py` is the runnable flagship example (build a
particle → topological features → `TOPCalculator` energy); `train_top.ipynb` covers fitting;
`multimet_go.ipynb` is the trimetallic (Pd/Au/Cu) Monte Carlo ordering-optimization workflow.

## Architecture

The standard workflow is a pipeline tied together by **string keys**:

```
build Nanoparticle  →  feature_classifier.compute_feature_vector(p)  →
energy_calculator.compute_energy(p)  →  optimization
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

- **`npl.optimization`** — chemical-ordering optimization for **two-element** particles via basin
  hopping. Subpackages: `basin_hopping`, `local_optimization`, and `genetic_algorithm`. Basin
  hopping / local optimization select the exchange operator via `model="TOP"` (guided,
  `GuidedExchangeOperator`) or `model="ACT"` (atomic coordination types, `ACTExchangeOperator`) —
  see `EXCHANGE_OPERATORS` in `local_optimization/local_optimization.py`. The guided exchange is
  driven by per-atom energies from the **TEC** local descriptor: for `model="TOP"`, fit/convert the
  global TOP (`TopologicalFeatureClassifier`) coefficients with
  `compute_coefficients_for_linear_topological_model` and load them onto a `TOPCalculator('TEC', …)`
  via `set_coefficients`; `model="ACT"` is already trained on the TEC descriptors directly. Guided
  exchange only exists for two elements, so **multi-metal (ETOP) ordering has no basin hopping —
  use the canonical Metropolis MC in `npl.monte_carlo` instead** (see `examples/multimet_go.ipynb`).

- **`npl.monte_carlo`** — Metropolis Monte Carlo over chemical ordering at fixed composition
  (canonical). `run_monte_carlo(temperature, max_steps, start_particle, energy_calculator,
  feature_classifier)` (in `monte_carlo_etop.py`) is the main entry point and the **only multi-metal
  ordering optimizer** — the basin-hopping / guided-exchange path in `npl.optimization` is restricted
  to two elements (`symbol_a`/`symbol_b`), so trimetallic+ ordering goes through MC. Uses random
  exchange operators with incremental feature updates (only the swapped atoms' neighborhood is
  recomputed). `monte_carlo.py`'s `mc_run` is an older adsorbate-MC variant.

- **`npl.visualize`**, **`npl.utils`** — plotting helpers and math/geometry utilities.

## Gotchas

- **`npl.monte_carlo` is canonical MC over chemical ordering** (fixed composition). It is distinct
  from the companion repo `mcpy` (https://github.com/farrisric/mcpy), which does *grand-canonical*
  MC with machine-learning interatomic potentials (variable particle number, MACE/alchemi). New
  GCMC / ML-potential work belongs in `mcpy`; chemical-ordering MC stays here. mcpy's calculators
  live at the additional working directory `/home/energystorage/projects/mcpy/mcpy/calculators`.
- `import npl` runs `setup_logging()`, which writes an `npl.log` file in the current working
  directory and logs to stdout at INFO. Several modules also call `logging.basicConfig` again.
- `TOPCalculator` reads its coefficients from the Python dict `top_parameters` in
  `npl/calculators/parameters.py` (keyed by stoichiometry string, e.g. `"Pt151Cu50"`).
- The version is single-sourced from `npl/__init__.py:__version__` (`pyproject.toml` reads it via
  `[tool.setuptools.dynamic]`). Bump it there.
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