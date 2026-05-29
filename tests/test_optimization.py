import numpy as np
import pytest
from npl.calculators import BayesianRRCalculator
from npl.core import Nanoparticle
from npl.calculators import (
    TOPCalculator,
    compute_coefficients_for_linear_topological_model,
)
from npl.descriptors import ExtendedTopologicalFeaturesClassifier
from npl.optimization import local_optimization
from npl.optimization.basin_hopping.basin_hopping import run_basin_hopping
from npl.optimization.local_optimization.etop_exchange_operator import ETOPExchangeOperator


def _build():
    symbols = ["Pt", "Cu"]
    calc = TOPCalculator(
        "ETOP", stoichiometry="Pt151Cu50",
        feature_classifier=ExtendedTopologicalFeaturesClassifier,
    )
    fc = calc.get_feature_classifier()
    particle = Nanoparticle()
    particle.truncated_octahedron(7, 2, {"Pt": 151, "Cu": 50})
    fc.compute_feature_vector(particle)
    env_energies, _ = compute_coefficients_for_linear_topological_model(
        calc.get_coefficients(), symbols, n_atoms=201
    )
    return particle, calc, env_energies


@pytest.mark.parametrize("model", ["TOP", "ACT"])
def test_local_optimization_runs_and_does_not_worsen(model):
    particle, calc, env_energies = _build()
    result_particle, accepted = local_optimization(
        particle, calc, env_energies, model=model
    )
    assert len(accepted) >= 1
    assert accepted[-1][0] <= accepted[0][0] + 1e-6


@pytest.mark.parametrize("model", ["TOP", "ACT"])
def test_basin_hopping_runs_for_both_models(model):
    particle, calc, env_energies = _build()
    best, lowest_energies, _ = run_basin_hopping(
        particle, calc, env_energies,
        n_hopping_attempts=2, n_hops=2, model=model,
    )
    assert best is not None
    assert len(lowest_energies) >= 1


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
