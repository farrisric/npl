from itertools import product

import numpy as np
import pytest
from npl.core import Nanoparticle
from npl.calculators import (
    TOPCalculator,
    BayesianRRCalculator,
    compute_coefficients_for_linear_topological_model,
)
from npl.descriptors import ExtendedTopologicalFeaturesClassifier
from npl.optimization import local_optimization
from npl.optimization.basin_hopping.basin_hopping import run_basin_hopping
from npl.optimization.local_optimization.local_optimization import (
    setup_local_optimization,
    update_atomic_features,
)
from npl.optimization.local_optimization.etop_exchange_operator import ETOPExchangeOperator


def _build():
    """Wire the guided two-element optimizer the way GO_tab_data.ipynb does:
    a ``TOPCalculator`` reading the per-atom ``TEC`` descriptor, with the global
    TOP coefficients converted into per-atom coefficients (``set_coefficients``)
    and the matching ``total_energies`` driving the guided exchange. The local
    optimization loop updates ``TEC`` incrementally, so the calculator MUST read
    ``TEC`` (not the global ``ETOP``) or the energy never changes."""
    np.random.seed(3)
    symbols = ["Pt", "Cu"]
    calc = TOPCalculator(
        "TEC", stoichiometry="Pt151Cu50",
        feature_classifier=ExtendedTopologicalFeaturesClassifier,
    )
    coefficients, total_energies = compute_coefficients_for_linear_topological_model(
        calc.get_coefficients(), symbols, n_atoms=201
    )
    calc.set_coefficients(coefficients)
    particle = Nanoparticle()
    particle.truncated_octahedron(7, 2, {"Pt": 151, "Cu": 50})
    return particle, calc, total_energies


@pytest.mark.parametrize("model", ["TOP", "ACT"])
def test_local_optimization_lowers_energy(model):
    particle, calc, total_energies = _build()
    result_particle, accepted = local_optimization(
        particle, calc, total_energies, model=model
    )
    start_energy, final_energy = accepted[0][0], accepted[-1][0]
    # a random ordering is never optimal: the guided loop must take real steps
    # and strictly lower the energy (the old ETOP wiring was a no-op, len == 2).
    assert len(accepted) > 2
    assert final_energy < start_energy


@pytest.mark.parametrize("model", ["TOP", "ACT"])
def test_basin_hopping_lowers_energy(model):
    particle, calc, total_energies = _build()
    best, lowest_energies, _ = run_basin_hopping(
        particle, calc, total_energies,
        n_hopping_attempts=2, n_hops=2, model=model,
    )
    assert best is not None
    assert lowest_energies[-1][0] < lowest_energies[0][0]


def _instrumented_descent(model):
    """Greedy guided descent recording, per accepted step, a tuple of
    ``(chosen_swap_energy, brute_force_min_over_all_pairs, actual_delta_E,
    swapped_atoms_are_neighbors)``. Mirrors ``local_optimization`` but exposes
    the exchange operator's predictions so they can be checked against the
    recomputed energy."""
    particle, calc, total_energies = _build()
    energy_key, lec, lfc, op = setup_local_optimization(
        particle, calc, total_energies, model=model
    )
    start_energy = particle.get_energy(energy_key)

    records = []
    while True:
        e1, e2 = op.symbol1_exchange_energies, op.symbol2_exchange_energies
        brute_min = min(e1[a] + e2[b] for a, b in product(e1, e2))

        i1, i2, _ = op.guided_exchange(particle)
        chosen = e1[i1] + e2[i2]
        neighbors = i2 in set(particle.neighbor_list[i1])

        particle, neighborhood = update_atomic_features(i1, i2, lec, lfc, particle)
        op.update(particle, neighborhood, [i1, i2])
        calc.compute_energy(particle)
        new_energy = particle.get_energy(energy_key)

        records.append((chosen, brute_min, new_energy - start_energy, neighbors))
        if new_energy < start_energy:
            start_energy = new_energy
        else:
            particle.swap_symbols([(i1, i2)])
            break
    return records


@pytest.mark.parametrize("model", ["TOP", "ACT"])
def test_guided_exchange_selects_lowest_predicted_swap(model):
    # The guided operator takes the head of each sorted per-atom list; since the
    # predicted swap energy is the separable sum e1[i1] + e2[i2], that must equal
    # the brute-force minimum over every symbol1 x symbol2 pair.
    records = _instrumented_descent(model)
    assert len(records) > 1
    for chosen, brute_min, _, _ in records:
        assert chosen == pytest.approx(brute_min, abs=1e-9)


def test_top_prediction_matches_actual_energy_for_nonneighbor_swaps():
    # For the TOP guided operator the predicted swap energy equals the true
    # change in energy exactly, as long as the two swapped atoms are not
    # neighbors (the per-atom prediction can't see their shared bond change).
    records = _instrumented_descent("TOP")
    nonneighbor = [(chosen, actual) for chosen, _, actual, nb in records if not nb]
    assert len(nonneighbor) > 0
    for chosen, actual in nonneighbor:
        assert chosen == pytest.approx(actual, abs=1e-6)


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
        # SortedKeyList order is insertion-order-dependent for tied gains, so
        # compare the sorted gain sequences rather than the index sequences.
        op_gains = [op.flip_gain[pair][idx] for idx in op.indices[pair]]
        fresh_gains = [fresh.flip_gain[pair][idx] for idx in fresh.indices[pair]]
        assert op_gains == pytest.approx(fresh_gains)


def test_etop_basin_hop_step_swaps_unlike_species():
    particle, calc, fc = _build_etop()
    op = ETOPExchangeOperator(calc.get_coefficients(), fc)
    op.bind_particle(particle)
    a, b = op.basin_hop_step(particle)
    assert op.atom_symbol[a] != op.atom_symbol[b]
    assert particle.get_symbol(a) == op.atom_symbol[b]


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


def test_etop_basin_hopping_runs():
    particle, calc, fc = _build_etop()
    best, lowest_energies, flip_list = run_basin_hopping(
        particle, calc, None, n_hopping_attempts=2, n_hops=2, model="ETOP")
    assert best is not None
    assert len(lowest_energies) >= 1
    assert lowest_energies[-1][0] <= lowest_energies[0][0] + 1e-6


# --------------------------------------------------------------------- GA
def _ga_setup(n_population=6):
    """A start population of random Pt151Cu50 orderings plus the wiring the GA
    needs. Same calculator as the guided-optimizer tests above."""
    from npl.descriptors import (
        NeighborCountingEnvironmentCalculator,
        TopologicalEnvironmentClassifier,
    )

    np.random.seed(5)
    _particle, calc, total_energies = _build()
    symbols = ["Pt", "Cu"]
    env_calc = NeighborCountingEnvironmentCalculator(symbols)
    classifier = TopologicalEnvironmentClassifier(env_calc, symbols)

    population = []
    for _ in range(n_population):
        p = Nanoparticle()
        p.truncated_octahedron(7, 2, {"Pt": 151, "Cu": 50})
        population.append(p)
    return population, calc, env_calc, classifier, total_energies


@pytest.mark.parametrize("model", ["ACT", "TOP"])
def test_ga_lowers_energy_and_keeps_composition(model):
    from npl.optimization.genetic_algorithm import run_single_particle_ga

    population, calc, env_calc, classifier, total_energies = _ga_setup()
    best_energies, best, evaluations = run_single_particle_ga(
        population, 3, calc, env_calc, classifier, total_energies, model=model,
    )
    assert best_energies[-1][0] < best_energies[0][0]
    assert evaluations > 0
    symbols = best.get_symbols()
    assert list(symbols).count("Pt") == 151
    assert list(symbols).count("Cu") == 50


def test_ga_defaults_to_the_act_operator():
    """'TOP' is deterministic in both its guided step and its kick, so every
    offspring collapses into the same limit cycle. The default must not be TOP."""
    import inspect

    from npl.optimization.genetic_algorithm import run_single_particle_ga

    assert inspect.signature(run_single_particle_ga).parameters['model'].default == 'ACT'


def test_ga_convergence_ignores_float_noise():
    """An improvement of 1e-12 is float noise. Counting it as progress kept the
    loop from ever converging."""
    from npl.optimization.genetic_algorithm import run_single_particle_ga

    population, calc, env_calc, classifier, total_energies = _ga_setup()
    _best_energies, _best, evaluations = run_single_particle_ga(
        population, 2, calc, env_calc, classifier, total_energies,
        improvement_tol=1e-6,
    )
    assert evaluations > 0     # terminated rather than spinning forever


def test_ga_uniform_fitness_when_population_is_degenerate():
    """All-equal energies must not zero every fitness: that used to break the
    generation loop and read as convergence."""
    from npl.optimization.genetic_algorithm import compute_fitness

    particle = Nanoparticle()
    particle.truncated_octahedron(7, 2, {"Pt": 151, "Cu": 50})
    particle.set_energy("TOP", -1.0)
    assert compute_fitness(particle, -1.0, -1.0, "TOP") > 0


def test_ordering_key_separates_degenerate_orderings():
    """Two different orderings can share an energy exactly; de-duplicating on
    energy would discard one of them."""
    from npl.optimization.genetic_algorithm import ordering_key

    a = Nanoparticle()
    a.truncated_octahedron(7, 2, {"Pt": 151, "Cu": 50})
    b = Nanoparticle()
    b.truncated_octahedron(7, 2, {"Pt": 151, "Cu": 50})
    indices_pt = a.get_indices_by_symbol("Pt")
    indices_cu = a.get_indices_by_symbol("Cu")
    b.swap_symbols([(int(indices_pt[0]), int(indices_cu[0]))])
    assert ordering_key(a) != ordering_key(b)


def test_cut_and_splice_preserves_composition_and_atom_count():
    from npl.optimization.genetic_algorithm import CutAndSpliceOperator

    np.random.seed(11)
    parents = []
    for _ in range(2):
        p = Nanoparticle()
        p.truncated_octahedron(7, 2, {"Pt": 151, "Cu": 50})
        parents.append(p)
    child = CutAndSpliceOperator().cut_and_splice(parents[0], parents[1])
    symbols = list(child.get_symbols())
    assert len(symbols) == 201
    assert symbols.count("Pt") == 151 and symbols.count("Cu") == 50
