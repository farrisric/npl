from npl.core import Nanoparticle
from npl.global_optimization.operations import GuidedExchangeOperator
from npl.global_optimization.operations import features
from npl.descriptors import Topologies
from npl.utils import compute_coefficients_for_linear_topological_model, get_combinations
from ase.cluster import Octahedron
from sortedcontainers import SortedKeyList
import random
from itertools import combinations_with_replacement
import numpy as np

atoms = Octahedron('Pt', 6, 1)

for i in range(20):
    atoms[i].symbol = 'Au'
    atoms[i+50].symbol = 'Ni'

system = Nanoparticle.from_atoms(atoms)

global_top = [-random.random() for _ in range(48)]
atomic_numbers = system.get_unique_atomic_numbers()

top = Topologies(atomic_numbers)
top.get_feature_vector(system)

coefficeints, feature_index_values = compute_coefficients_for_linear_topological_model(global_top, atomic_numbers)
exchange = GuidedExchangeOperator(system, coefficeints, feature_index_values)

sorted_energies_dict = {Z : SortedKeyList(key=lambda x: min(exchange.exchange_energies[x])) for Z in exchange.atomic_numbers}

for i in exchange.sorted_indices:
    Z = system[i].number
    sorted_energies_dict[Z].add(i)

atomic_numbers = {i : x for i, x in enumerate(exchange.atomic_numbers)}

exchange_matrix = np.empty((3,3))

for f in range(2):
    for Z_index, Z in atomic_numbers.items():
        i = sorted_energies_dict[Z][f]
        print(i)
        flip = exchange.exchange_energies[i]
        exchange_matrix[Z_index] = flip

    combinations = get_combinations(range(3),2)
    combinations = [x for x in combinations if x[1]!=x[0]]

    for combination in combinations:
        i,j = combination 
        exchange_ij = exchange_matrix[i][j] + exchange_matrix[j][i]
        print(combination ,exchange_ij)

