"""The Rust ACT fast path must agree with the Python it replaces.

Three things are checked, in rising order of how much they would cost to get
wrong:

1. the energy model - the Rust sum over ``coef[feature]`` equals the Python
   dot product of the ridge coefficients with the TEC feature vector;
2. the incremental exchange-energy bookkeeping - after a full search the
   maintained exchange energies still equal a from-scratch recompute. This is
   the part of the port that is easy to get subtly wrong, and it is checked
   independently of which atom a tie picks;
3. the search itself - the deterministic descent reaches the same energy as
   the Python descent, and basin hopping never returns a particle whose
   reported energy disagrees with the Python calculator.
"""

import numpy as np
import pytest

from npl.optimization import act_fast
from npl.optimization.basin_hopping.basin_hopping import run_basin_hopping
from npl.optimization.local_optimization.local_optimization import local_optimization

from test_optimization import _build  # noqa: E402  (shared fixture builder)


pytestmark = pytest.mark.skipif(not act_fast.available(),
                                reason='npl_act_rs is not installed')


def _setup():
    """A Pt/Cu particle, a fitted linear TEC calculator and its environment
    energies - the same fixture the existing optimization tests use."""
    particle, calculator, environment_energies = _build()
    return particle, calculator, environment_energies


def test_energy_matches_the_python_calculator():
    particle, calculator, _ = _setup()
    from npl.descriptors import (NeighborCountingEnvironmentCalculator,
                                 TopologicalEnvironmentClassifier)

    symbols = particle.get_all_symbols()
    env_calculator = NeighborCountingEnvironmentCalculator(symbols)
    classifier = TopologicalEnvironmentClassifier(env_calculator, symbols)
    env_calculator.compute_local_environments(particle)
    classifier.compute_feature_vector(particle)

    expected = calculator.compute_energy(particle)

    import npl_act_rs
    offsets, neighbours, codes, _elements = act_fast.pack(particle)
    got = npl_act_rs.act_energy(
        offsets, neighbours, codes,
        [float(c) for c in calculator.get_coefficients()])

    assert got == pytest.approx(expected, abs=1e-9)


def test_packing_round_trips():
    particle, _, _ = _setup()
    before = np.array(particle.atoms.atoms.numbers)
    offsets, _neighbours, codes, elements = act_fast.pack(particle)

    assert len(offsets) == len(before) + 1
    assert offsets[0] == 0
    act_fast.unpack(particle, codes, elements)
    assert np.array_equal(particle.atoms.atoms.numbers, before)


def test_local_optimization_reaches_the_same_energy():
    """The descent is deterministic in both implementations, so they must land
    at the same energy. The ordering itself may differ where two atoms carry
    equal exchange energies and the tie is broken differently."""
    _particle, calculator, environment_energies = _setup()

    rust_particle, rust_accepted = local_optimization(
        _setup()[0], calculator, environment_energies, model='ACT')
    python_particle, python_accepted = local_optimization(
        _setup()[0], calculator, environment_energies, model='ACT',
        force_python=True)

    rust_energy = rust_accepted[-1][0]
    python_energy = python_accepted[-1][0]
    assert rust_energy == pytest.approx(python_energy, abs=1e-6)

    # and the reported energy is really the energy of the returned ordering
    assert calculator.compute_energy(rust_particle) == pytest.approx(
        rust_energy, abs=1e-6)
    assert python_particle is not None


def test_basin_hopping_energy_is_consistent_and_no_worse():
    """Basin hopping is stochastic and the two RNG streams differ, so this
    checks what must hold for any correct implementation: the reported energy
    belongs to the returned ordering, the run improves on its start, and the
    Rust result is in the same league as the Python one."""
    particle, calculator, environment_energies = _setup()
    start_energy = None

    rust_best, rust_lowest, rust_flips = run_basin_hopping(
        _setup()[0], calculator, environment_energies,
        n_hopping_attempts=5, n_hops=3, model='ACT', seed=12345)
    python_best, python_lowest, _ = run_basin_hopping(
        _setup()[0], calculator, environment_energies,
        n_hopping_attempts=5, n_hops=3, model='ACT', force_python=True)

    start_energy = rust_lowest[0][0]
    assert rust_lowest[-1][0] <= start_energy + 1e-9
    assert calculator.compute_energy(rust_best) == pytest.approx(
        rust_lowest[-1][0], abs=1e-6)
    assert len(rust_flips) > 0
    # every flip entry is (predicted flip, realised dE, [i, j])
    assert all(len(pair) == 2 for _f, _d, pair in rust_flips)

    # same search, same landscape: within a few meV of each other
    assert rust_lowest[-1][0] == pytest.approx(python_lowest[-1][0], abs=0.05)


def test_bookkeeping_stays_exact_through_a_search():
    """After a full basin-hopping run the incrementally updated exchange
    energies must still equal a from-scratch recompute."""
    particle, calculator, environment_energies = _setup()
    _result, error = act_fast.basin_hopping(
        particle, calculator, environment_energies,
        n_hopping_attempts=5, n_hops=3,
        local_env_calculator=_env_calculator(particle),
        local_feature_classifier=_classifier(particle),
        seed=99)
    assert error < 1e-9


def _env_calculator(particle):
    from npl.descriptors import NeighborCountingEnvironmentCalculator
    return NeighborCountingEnvironmentCalculator(particle.get_all_symbols())


def _classifier(particle):
    from npl.descriptors import TopologicalEnvironmentClassifier
    return TopologicalEnvironmentClassifier(_env_calculator(particle),
                                            particle.get_all_symbols())


def test_falls_back_when_the_model_is_not_act():
    particle, calculator, environment_energies = _setup()
    assert not act_fast.usable(particle, calculator, environment_energies, 'TOP')
