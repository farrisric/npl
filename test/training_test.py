from random import shuffle

from ase.cluster import Octahedron
from ase.calculators.emt import EMT
from ase.visualize import view

from sklearn.linear_model import Ridge

import matplotlib.pyplot as plt

from npl.calculators import TrainedCalculator
from npl.descriptors import Topologies
from npl.core import Nanoparticle as nano
from npl.utils import compute_coefficients_for_linear_topological_model
from npl.global_optimization.operations import GuidedExchangeOperator

symbols = ['Au']*150+['Pt']*100+['Pd']*155
system = Octahedron('Pt', 9, 3)
system.calc = EMT()
training_set = []

for _ in range(20):
    shuffle(symbols)
    system.symbols = symbols
    system.get_potential_energy()
    p = nano.from_atoms(system)
    top = Topologies(p.get_unique_atomic_numbers())
    top.get_feature_vector(p)
    training_set.append(p)

c = TrainedCalculator('TOP', Ridge)
c.fit(training_set)
coef = c.get_coefficients()

coefficients, feature_index_values = compute_coefficients_for_linear_topological_model(coef, training_set[1].get_unique_atomic_numbers())
exchange = GuidedExchangeOperator(p, coefficients, feature_index_values)

start_energy = c.calculate_total(p)

exchanges = []
energies = []
for _ in range(10):
    best_swap = exchange.perform_operation(p)
    top.get_feature_vector(p)
    swap = c.calculate_total(p) - start_energy
    start_energy = c.calculate_total(p)
    print(best_swap, swap)
    exchanges.append(best_swap)
    energies.append(swap)


plt.scatter(energies, exchanges)
plt.show()