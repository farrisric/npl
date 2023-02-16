from npl.core import Nanoparticle
from npl.global_optimization.operations import GuidedExchangeOperator
from npl.global_optimization.operations import features
from npl.descriptors import Topologies
from npl.utils import compute_coefficients_for_linear_topological_model

from ase.cluster import Octahedron
from sortedcontainers import SortedKeyList
import random
from itertools import combinations_with_replacement
import numpy as np

atoms = Octahedron('Pt', 6, 1)
atoms[10].symbol = 'Au'
atoms[12].symbol = 'Ni'

system = Nanoparticle.from_atoms(atoms)

global_top = [-random.random() for _ in range(48)]
symbols = ['Au','Pt','Ni']

top = Topologies(symbols)

coefficeints, feature_index_values = compute_coefficients_for_linear_topological_model(global_top, symbols)
exchange = GuidedExchangeOperator(system, coefficeints, feature_index_values)

top.get_feature_vector(system)

a = exchange.get_atom_feature(system, system[10])
print(a)