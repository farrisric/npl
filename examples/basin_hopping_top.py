"""Train a TOP model, then optimize chemical ordering with basin hopping.

Fits a linear TOP energy model on EMT reference data and then optimizes the
chemical ordering of a bimetallic Pt/Cu truncated octahedron with guided basin
hopping (``model="TOP"``).

The TOP model is fit on the *global* topological descriptor
(``testTopologicalFeatureClassifier`` -- the bond-count + coordination layout
the linear model expects). For basin hopping the fitted coefficients are
converted into per-atom environment energies with
``compute_coefficients_for_linear_topological_model``, and the calculator is
switched to read the per-atom ``TEC`` descriptor that the guided exchange
updates incrementally.

Guided exchange exists only for two elements; for multi-metal (ETOP) ordering
use the Metropolis Monte Carlo workflow in ``examples/multimet_go.ipynb``.
"""
from npl.core import Nanoparticle
from npl.calculators import (
    EMTCalculator,
    BayesianRRCalculator,
    compute_coefficients_for_linear_topological_model,
)
from npl.descriptors.global_feature_classifier import testTopologicalFeatureClassifier
from npl.optimization import run_basin_hopping

STOICHIOMETRY = {"Pt": 40, "Cu": 39}  # 79-atom truncated octahedron (height=5, trunc=1)
N_ATOMS = sum(STOICHIOMETRY.values())
SYMBOLS = list(STOICHIOMETRY)


def make_particle():
    particle = Nanoparticle()
    particle.truncated_octahedron(5, 1, STOICHIOMETRY)
    return particle


def build_training_set(n_particles):
    emt = EMTCalculator(fmax=0.4)
    training_set = []
    for _ in range(n_particles):
        particle = make_particle()
        emt.compute_energy(particle)
        training_set.append(particle)
    return training_set


def train_top_model(training_set):
    """Fit the linear TOP model on the global topological descriptor and return
    a calculator that evaluates the per-atom (TEC) energy together with the
    per-environment energies that drive the guided exchange."""
    classifier = testTopologicalFeatureClassifier(SYMBOLS)
    for particle in training_set:
        classifier.compute_feature_vector(particle)

    calculator = BayesianRRCalculator(classifier.get_feature_key())
    calculator.fit(training_set, "EMT")

    coefficients, total_energies = compute_coefficients_for_linear_topological_model(
        calculator.get_coefficients(), SYMBOLS, n_atoms=N_ATOMS
    )
    calculator.set_coefficients(coefficients)
    calculator.set_feature_key("TEC")
    return calculator, total_energies


def main():
    training_set = build_training_set(20)
    calculator, total_energies = train_top_model(training_set)

    particle = make_particle()
    best_particle, lowest_energies, _ = run_basin_hopping(
        particle, calculator, total_energies,
        n_hopping_attempts=5, n_hops=2, model="TOP",
    )

    start_energy, best_energy = lowest_energies[0][0], lowest_energies[-1][0]
    print(f"TOP basin hopping: {start_energy:.4f} -> {best_energy:.4f}")
    return best_particle


if __name__ == "__main__":
    main()
