"""Monte Carlo chemical-ordering optimization.

`npl.monte_carlo.run_monte_carlo` is the multi-metal ordering optimizer: the
basin-hopping / guided-exchange path is restricted to two elements, so Metropolis
MC over a fitted energy model is the route for trimetallic (and beyond) particles.
"""
import numpy as np

from npl.core import Nanoparticle
from npl.calculators import EMTCalculator, BayesianRRCalculator
from npl.descriptors import ExtendedTopologicalFeaturesClassifier
from npl.monte_carlo import run_monte_carlo


def _fit_brr(stoichiometry, n_train=15, height=4, trunc=1):
    classifier = ExtendedTopologicalFeaturesClassifier(list(stoichiometry.keys()))
    emt = EMTCalculator(fmax=0.2)
    training_set = []
    for _ in range(n_train):
        p = Nanoparticle()
        p.truncated_octahedron(height, trunc, stoichiometry)
        emt.compute_energy(p)
        classifier.compute_feature_vector(p)
        training_set.append(p)
    calc = BayesianRRCalculator(classifier.get_feature_key())
    calc.fit(training_set, "EMT")
    return calc, classifier, height, trunc


def test_run_monte_carlo_trimetallic_runs_and_conserves_composition():
    stoichiometry = {"Pd": 0.34, "Au": 0.33, "Cu": 0.33}
    calc, classifier, height, trunc = _fit_brr(stoichiometry)

    start = Nanoparticle()
    start.truncated_octahedron(height, trunc, stoichiometry)
    start_stoich = dict(start.get_stoichiometry())

    best, accepted = run_monte_carlo(300, 100, start, calc, classifier)

    assert best is not None
    assert len(accepted) >= 1
    assert np.isfinite(accepted[-1][0])
    # canonical MC: composition is conserved by the exchange moves
    assert dict(best.get_stoichiometry()) == start_stoich
