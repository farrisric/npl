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
    # best_particle must be the lowest-energy configuration actually sampled
    min_accepted = min(e for e, _ in accepted)
    assert abs(best.get_energy(calc.get_energy_key()) - min_accepted) < 1e-9


def test_run_monte_carlo_is_deterministic_given_seed():
    # Seeding *before* building the particle fixes the random initial ordering as
    # well as the Metropolis walk, so the whole run is reproducible. Guards against
    # changes that perturb the RNG stream or the feature-update order (e.g. the
    # single-pass symbol-index lookup in the random exchange operator).
    def one_run():
        np.random.seed(123)
        classifier = ExtendedTopologicalFeaturesClassifier(["Au", "Cu", "Pt"])
        p = Nanoparticle()
        p.truncated_octahedron(4, 1, {"Pt": 0.34, "Au": 0.33, "Cu": 0.33})
        classifier.compute_feature_vector(p)
        coeffs = np.random.default_rng(0).standard_normal(classifier.n_features)
        calc = BayesianRRCalculator("ETOP")
        calc.set_coefficients(coeffs)
        best, accepted = run_monte_carlo(300, 100, p, calc, classifier)
        return accepted, [int(x) for x in best.atoms.atoms.numbers], best.get_energy("BRR")

    accepted1, numbers1, energy1 = one_run()
    accepted2, numbers2, energy2 = one_run()
    assert accepted1 == accepted2
    assert numbers1 == numbers2
    assert energy1 == energy2
