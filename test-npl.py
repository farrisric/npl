from npl.descriptors import ExtendedTopologicalFeaturesClassifier
from npl.calculators import TOPCalculator
from npl.core import Nanoparticle
from npl.monte_carlo import run_monte_carlo
from npl.visualize import plot_parted_particle

energy_calculator = TOPCalculator('ETOP', stoichiometry='Pt151Cu50',
                                  feature_classifier=ExtendedTopologicalFeaturesClassifier)

feature_classifier = energy_calculator.get_feature_classifier()

temperature = 250
max_steps = 10000

start_particle = Nanoparticle()
start_particle.truncated_octahedron(7, 2, {'Pt': 151, 'Cu': 50})
best_particle, accepted_energies = run_monte_carlo(temperature,
                                                   max_steps,
                                                   start_particle,
                                                   energy_calculator,
                                                   feature_classifier)

plot_parted_particle(best_particle)
