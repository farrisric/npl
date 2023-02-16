from npl.core import Nanoparticle
from npl.global_optimization.operations import GuidedExchangeOperator
from npl.utils import compute_coefficients_for_linear_topological_model
from ase.cluster import Octahedron
import random

atoms = Octahedron('Pt', 6, 1)
atoms[10].symbol = 'Au'
atoms[12].symbol = 'Ni'

system = Nanoparticle.from_atoms(atoms)

global_top = [-random.random() for _ in range(48)]
symbols = ['Au','Pt','Ni']

coefficeints, feature_index_values = compute_coefficients_for_linear_topological_model(global_top, symbols)
exchange = GuidedExchangeOperator(system, coefficeints)

