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


def test_local_optimization_runs_and_does_not_worsen():
    particle, calc, env_energies = _build()
    result_particle, accepted = local_optimization(particle, calc, env_energies)
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
