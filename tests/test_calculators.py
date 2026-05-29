import numpy as np
from npl.core import Nanoparticle
from npl.calculators import TOPCalculator, EMTCalculator, BayesianRRCalculator
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


def test_brr_fit_float_validation_split_uses_majority_for_training():
    """A float validation_set must reserve that *fraction* for validation and
    train on the rest. Regression test: the split was inverted, training on the
    small fraction (and on an empty set for <=9 particles)."""
    stoichiometry = {"Pd": 0.34, "Au": 0.33, "Cu": 0.33}
    etop = ExtendedTopologicalFeaturesClassifier(list(stoichiometry.keys()))

    training_set = []
    for i in range(8):
        p = Nanoparticle()
        p.truncated_octahedron(5, 1, stoichiometry)
        etop.compute_feature_vector(p)
        p.set_energy("EMT", float(i))
        training_set.append(p)

    calc = BayesianRRCalculator(etop.get_feature_key())

    captured = {}
    real_fit = calc.ridge.fit

    def spy_fit(X, y):
        captured["n_train"] = len(X)
        return real_fit(X, y)

    calc.ridge.fit = spy_fit
    calc.fit(training_set, "EMT", validation_set=0.1)

    # 8 particles, 10% validation -> 7 train / 1 validation (not 0 train / 8 val)
    assert captured["n_train"] == 7
