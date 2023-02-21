import random
from itertools import combinations_with_replacement

import numpy as np
from ase.cluster import Octahedron
from sortedcontainers import SortedKeyList

from npl.core import Nanoparticle
from npl.descriptors import Topologies
from npl.global_optimization.operations import GuidedExchangeOperator, features
from npl.utils import (compute_coefficients_for_linear_topological_model,
                       get_combinations)

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

exchange = exchange.find_best_swap_pair()
print(exchange)

# index_atomic_numbers = {x : i for i, x in enumerate(exchange.atomic_numbers)}

# sorted_energies_dict = {Z : {} for Z in index_atomic_numbers.values()}

# for Z_i in sorted_energies_dict:
#     for Z_j in sorted_energies_dict:
#         if Z_j == Z_i:
#             continue
#         sorted_energies_dict[Z_i][Z_j] = SortedKeyList(key=lambda x: exchange.exchange_energies[x][Z_j])

# for atom in system:
#     Z_i = index_atomic_numbers[atom.number]
#     for Z_j in index_atomic_numbers.values():
#         if Z_i == Z_j:
#             continue
#         sorted_energies_dict[Z_i][Z_j].add(atom.index)

# combinations = get_combinations(range(3),2)
# combinations = [x for x in combinations if x[1]!=x[0]]

# best_flip = 0
# for combination in combinations:
#     Z_i, Z_j = combination
#     a = sorted_energies_dict[Z_i][Z_j][0]
#     b = sorted_energies_dict[Z_j][Z_i][0]
#     flip = exchange.exchange_energies[a][Z_j] + exchange.exchange_energies[b][Z_i]
#     if flip < best_flip:
#         exchange_indices = (a, b)
#         best_flip = flip

# print(best_flip, exchange_indices)
# print(system[exchange_indices[0]])
# print(system[exchange_indices[1]])
