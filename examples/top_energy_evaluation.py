"""Flagship example: evaluate the topological (TOP) energy of a bimetallic
nanoparticle using pretrained coefficients.

Builds a 201-atom truncated-octahedral Pt151Cu50 particle, computes its
extended topological feature vector, and evaluates the energy with a
TOPCalculator loaded from pretrained coefficients shipped with npl.
"""
from npl.core import Nanoparticle
from npl.calculators import TOPCalculator
from npl.descriptors import ExtendedTopologicalFeaturesClassifier


def main():
    calculator = TOPCalculator(
        "ETOP",
        stoichiometry="Pt151Cu50",
        feature_classifier=ExtendedTopologicalFeaturesClassifier,
    )
    feature_classifier = calculator.get_feature_classifier()

    particle = Nanoparticle()
    particle.truncated_octahedron(7, 2, {"Pt": 151, "Cu": 50})

    feature_classifier.compute_feature_vector(particle)
    energy = calculator.compute_energy(particle)

    print(f"TOP energy of Pt151Cu50: {float(energy):.4f}")
    return float(energy)


if __name__ == "__main__":
    main()
