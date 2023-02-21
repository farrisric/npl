import random
from itertools import combinations_with_replacement

import numpy as np
from ase.cluster import Octahedron
from ase.visualize import view
from sortedcontainers import SortedKeyList

from npl.core import Nanoparticle
from npl.descriptors import Topologies
from npl.global_optimization.operations import GuidedExchangeOperator, features
from npl.utils import (compute_coefficients_for_linear_topological_model,
                       get_combinations)

atoms = Octahedron('Pt', 9, 3)

for i in range(100):
    atoms[i].symbol = 'Au'
    atoms[i+50].symbol = 'Ni'

system = Nanoparticle.from_atoms(atoms)
view(system)

global_top = [-random.random() for _ in range(48)]
atomic_numbers = system.get_unique_atomic_numbers()

top = Topologies(atomic_numbers)
top.get_feature_vector(system)

coefficeints, feature_index_values = compute_coefficients_for_linear_topological_model(global_top, atomic_numbers)
exchange = GuidedExchangeOperator(system, coefficeints, feature_index_values)

for _ in range(1000):
    exchange.perform_operation(system)
view(system)