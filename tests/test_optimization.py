from itertools import product

import numpy as np
import pytest
from npl.core import Nanoparticle
from npl.calculators import (
    TOPCalculator,
    compute_coefficients_for_linear_topological_model,
)
from npl.descriptors import ExtendedTopologicalFeaturesClassifier
from npl.optimization import local_optimization
from npl.optimization.basin_hopping.basin_hopping import run_basin_hopping
from npl.optimization.local_optimization.local_optimization import (
    setup_local_optimization,
    update_atomic_features,
)


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

        i1, i2 = op.guided_exchange(particle)
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
