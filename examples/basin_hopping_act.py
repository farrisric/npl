"""Train an ACT model, then optimize chemical ordering with basin hopping.

Fits a linear ACT energy model on EMT reference data and then optimizes the
chemical ordering of a bimetallic Pt/Cu truncated octahedron with guided basin
hopping using the atomic-coordination-types exchange operator
(``model="ACT"``, ``ACTExchangeOperator``).

Unlike TOP (fit on the global descriptor and then converted), the ACT model is
fit *directly* on the per-atom topological environment descriptor (``TEC``).
The fitted linear coefficients ARE the per-environment energies, so the same
array is used both as the calculator coefficients and as the energies driving
the ACT exchange -- there is no conversion step.

Guided exchange exists only for two elements; for multi-metal (ETOP) ordering
use the Metropolis Monte Carlo workflow in ``examples/multimet_go.ipynb``.
"""
from npl.core import Nanoparticle
from npl.calculators import EMTCalculator, BayesianRRCalculator
from npl.descriptors import (
    NeighborCountingEnvironmentCalculator,
    TopologicalEnvironmentClassifier,
)
from npl.optimization import run_basin_hopping

STOICHIOMETRY = {"Pt": 40, "Cu": 39}  # 79-atom truncated octahedron (height=5, trunc=1)
SYMBOLS = list(STOICHIOMETRY)


def make_particle():
    particle = Nanoparticle()
    particle.truncated_octahedron(5, 1, STOICHIOMETRY)
    return particle


def build_training_set(n_particles, env_calculator, classifier):
    emt = EMTCalculator(fmax=0.4)
    training_set = []
    for _ in range(n_particles):
        particle = make_particle()
        emt.compute_energy(particle)
        env_calculator.compute_local_environments(particle)
        classifier.compute_feature_vector(particle)
        training_set.append(particle)
    return training_set


def train_act_model(training_set, classifier):
    """Fit the linear model directly on the per-atom TEC descriptor; the fitted
    coefficients are the per-environment energies used by the ACT exchange."""
    calculator = BayesianRRCalculator(classifier.get_feature_key())
    calculator.fit(training_set, "EMT")
    environment_energies = calculator.get_coefficients()
    return calculator, environment_energies


def main():
    env_calculator = NeighborCountingEnvironmentCalculator(SYMBOLS)
    classifier = TopologicalEnvironmentClassifier(env_calculator, SYMBOLS)

    training_set = build_training_set(20, env_calculator, classifier)
    calculator, environment_energies = train_act_model(training_set, classifier)

    particle = make_particle()
    best_particle, lowest_energies, _ = run_basin_hopping(
        particle, calculator, environment_energies,
        n_hopping_attempts=5, n_hops=2, model="ACT",
    )

    start_energy, best_energy = lowest_energies[0][0], lowest_energies[-1][0]
    print(f"ACT basin hopping: {start_energy:.4f} -> {best_energy:.4f}")
    return best_particle


if __name__ == "__main__":
    main()
