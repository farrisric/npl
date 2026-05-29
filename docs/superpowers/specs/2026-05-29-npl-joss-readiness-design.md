# NPL JOSS-Readiness — Design

**Date:** 2026-05-29
**Author:** Riccardo Farris (with Claude Code)
**Status:** Draft for review

## Goal

Bring NanoParticleLibrary (`npl`) to a state where it can be submitted to the
Journal of Open Source Software (JOSS): the library installs and runs on a
current scientific-Python stack, dead code and broken exports are gone, an
automated test suite runs in CI, packaging is modern, and a JOSS paper
(`paper.md` + `paper.bib`) is drafted.

mcpy (the sibling repo at https://github.com/farrisric/mcpy) is the quality
bar to mirror: it already has `pyproject.toml`, a real `tests/` suite,
`CONTRIBUTING.md`, `CLAUDE.md`, and a JOSS `paper/`. After this work, npl
should match that bar, and the npl paper should *reference* mcpy rather than
duplicate Monte Carlo functionality.

## Scope decisions (confirmed with user)

- **Monte Carlo:** remove `npl/monte_carlo` from npl entirely; it lives in mcpy
  now. Migrate examples/docs off it.
- **Tests:** build a full pytest suite and re-enable CI.
- **Paper:** draft the full JOSS `paper.md` + `paper.bib` this round.
- **Compatibility:** target the current stack only (numpy 2.x, scipy 1.17,
  ASE 3.28, scikit-learn 1.8, Python 3.10–3.13). No back-compat shims for old
  deps.
- **Packaging:** migrate from `setup.py` to PEP 621 `pyproject.toml` (mirror
  mcpy); single-source the version.
- **Flagship example:** replace the README's MC example with a *TOP energy
  evaluation* example — build a particle, compute topological features,
  evaluate energy with a pretrained `TOPCalculator`.

## Current state (audited 2026-05-29)

The library is **non-functional on a current stack**. On Python 3.13 /
numpy 2.4 / scipy 1.17 / ASE 3.28 / sklearn 1.8:

- `import npl.descriptors` (and transitively `calculators`, `optimization`,
  `monte_carlo`) fails: `npl/descriptors/local_environment_calculator.py:3`
  does `from scipy.special import sph_harm`, removed in modern SciPy
  (now `sph_harm_y`, with swapped argument order).
- `npl/core/nanoparticle.py:44` uses `np.int`, removed in NumPy 2.
- `npl/calculators/energy_calculator.py:104` uses `atoms.set_calculator(EMT())`,
  removed in ASE.

Only `npl.core`, `npl.visualize`, `npl.utils` import cleanly.

Pre-existing brokenness independent of deps:

- Broken `__all__` exports:
  - `npl/monte_carlo/__init__.py` exports `BaseEnsemble`, `CanonicalEnsemble`
    (do not exist).
  - `npl/descriptors/__init__.py` exports `LayererTopologicalDescriptors`
    (class exists in `global_feature_classifier.py` but is never imported into
    the package namespace).
  - `npl/optimization/__init__.py` has a missing comma:
    `"local_optimization" "run_basin_hopping"` silently concatenates into one
    string.
- Version mismatch: `setup.py` says `1.0.4`, `npl/__init__.py` says `1.0.0`.
- `npl/monte_carlo` has **zero internal importers** within npl (only
  examples/docs reference it), so removal is clean code-wise.
- `data/params.json` is unused — `TOPCalculator` reads coefficients from the
  Python dict `npl/calculators/parameters.py:top_parameters`.
- Dead-code suspects to confirm during implementation: the whole
  `npl/monte_carlo/` tree (to be removed), `garbage_exchange_operator.py`,
  `compute_coefficients_for_shape_optimization`, commented-out code blocks in
  `energy_calculator.py`, and the misspelled `genetic_algoritm` directory.

## Approach

Modernize-first, in verified vertical slices. Each phase has an exit criterion
that must pass before the next begins. (Tests-first was rejected: the code does
not import, so there is no running behavior to characterize first.)

## Phases

### Phase 1 — Restore importability on the current stack
- `sph_harm` → `sph_harm_y` in `local_environment_calculator.py`, correcting
  the argument order/convention (verify the spherical-harmonic values match the
  previous convention via a numerical check).
- `np.int` → `int` in `nanoparticle.py`.
- `atoms.set_calculator(EMT())` → `atoms.calc = EMT()` in `energy_calculator.py`.
- Sweep for any other removed/renamed APIs across numpy/scipy/ASE/sklearn.
- **Exit criterion:** `python -c "import npl; import npl.core, npl.descriptors,
  npl.calculators, npl.optimization, npl.visualize, npl.utils"` succeeds with no
  errors.

### Phase 2 — Remove `npl.monte_carlo`; redraw the boundary
- Delete `npl/monte_carlo/`.
- Fix broken `__all__` in `descriptors/__init__.py` and
  `optimization/__init__.py` (the comma bug and the missing
  `LayererTopologicalDescriptors` import — either import it or drop it from
  `__all__`, decided during implementation based on whether it is real/used).
- Migrate examples/docs that reference Monte Carlo:
  - Non-MC workflows rewritten against `npl.optimization`.
  - Genuine MC content replaced with a pointer to mcpy.
  - Affected: `README.md`, `examples/{global_optimization, multimet_go,
    loading_pretrain_top, adsorbate_optimization}.ipynb`, and the corresponding
    `docs/` and `docs/examples/` `.rst` files, plus `docs/npl.monte_carlo*.rst`.
- Replace the README flagship example with the TOP energy-evaluation example.
- **Exit criterion:** no remaining references to `npl.monte_carlo` in code,
  examples, or docs; the new flagship example runs end-to-end.

### Phase 3 — Dead-code sweep
- Remove commented-out code, unused functions, and broken/duplicate exports.
- Confirm-then-remove suspects (`garbage_exchange_operator.py`,
  `compute_coefficients_for_shape_optimization`, etc.); keep anything still
  referenced by examples, docs, or the public API.
- Remove the unused `data/params.json` (or wire it in) — decide during
  implementation; default is removal since `parameters.py` is the source of
  truth.
- Rename `genetic_algoritm` → `genetic_algorithm`, preserving the public import
  path used elsewhere (update all references).
- **Exit criterion:** `flake8` clean (config already present); no unused
  imports/exports; every removal traceable to this spec.

### Phase 4 — Test suite + CI
- Build a `pytest` suite under `test/` (rename to `tests/` to match mcpy) with
  unit coverage for:
  - `core`: particle construction, neighbor list, stoichiometry/ordering,
    save/load round-trip.
  - `descriptors`: feature-vector shape/values for the topological classifiers;
    local environment calculator (incl. the `sph_harm_y` fix).
  - `calculators`: `EMTCalculator`, `BayesianRRCalculator`, `TOPCalculator`
    (pretrained-coefficient load + energy evaluation).
  - `optimization`: a short basin-hopping / local-optimization run completes and
    lowers (or does not raise on) energy.
- Re-enable pytest in `.github/workflows/test.yml` on a Python matrix
  (3.10–3.13).
- **Exit criterion:** `pytest` passes locally; CI workflow runs pytest green.

### Phase 5 — Packaging + metadata
- Replace `setup.py` with PEP 621 `pyproject.toml` (mirror mcpy structure).
- Single-source the version (one canonical location; `__init__.__version__`
  derived or matched). Pick the target version during implementation
  (proposal: `1.1.0`, reflecting the MC removal + modernization).
- Align dependency floors with the current stack.
- Update `.github/workflows/publish.yml` to `python -m build`.
- **Exit criterion:** `python -m build` produces sdist+wheel; `pip install` of
  the wheel into a clean env imports `npl` and reports a single consistent
  version.

### Phase 6 — JOSS paper + repo essentials
- `paper/paper.md` (JOSS front-matter, Summary, Statement of Need, brief
  feature overview, acknowledgements) and `paper/paper.bib` containing the three
  references from the README plus a cross-citation to mcpy.
- `CONTRIBUTING.md` (mirror mcpy) and a `CITATION.cff` correctness check.
- Update `README.md` badges/links as needed.
- **Exit criterion:** `paper.md` has valid JOSS YAML front-matter and all
  `paper.bib` keys resolve; repo has CONTRIBUTING + accurate CITATION.cff.

## Success criteria (JOSS checklist mapped)

- Installs cleanly from a fresh environment; all public submodules import.
- Automated tests exist and run in CI.
- Example usage works and is documented.
- Statement of need + summary present in `paper.md`.
- Community guidelines (CONTRIBUTING) present.
- Versioning and packaging are consistent and modern.

## Risks / watch-items

- **`sph_harm_y` argument convention:** the comment at
  `local_environment_calculator.py:40` already notes a phi/theta convention
  quirk; the migration must preserve numerical results, verified by a test.
- **Example migration scope:** several notebooks lean on MC; rewriting them on
  `npl.optimization` may surface further breakage to fix in Phase 1/3.
- **Pretrained-model availability for the flagship example:** confirm a
  stoichiometry in `top_parameters` (e.g. a Pt/Au or Pt/Cu system) loads and
  evaluates so the README example is runnable and CI-testable.

## Non-goals

- No new scientific features or algorithms.
- No changes to mcpy.
- No re-derivation of the topological energy model or its coefficients.
