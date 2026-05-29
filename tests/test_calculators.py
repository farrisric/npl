import numpy as np
from npl.core import Nanoparticle
from npl.calculators import TOPCalculator, EMTCalculator
from npl.descriptors import ExtendedTopologicalFeaturesClassifier


def test_top_calculator_loads_and_evaluates(ptcu_particle):
    calc = TOPCalculator(
        "ETOP", stoichiometry="Pt151Cu50",
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
