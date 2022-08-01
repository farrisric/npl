import numba
from numba.experimental import jitclass
import numpy as np

from npl.descriptors import Descriptor
from npl.core import Nanoparticle

class EnvironmentalTopologies(Descriptor):
    
    def __init__(self, n_atoms):
        self.n_atoms = n_atoms
        self.local_environments = np.empty((n_atoms,2))
        self.atom_features = np.empty(n_atoms)

        self.coordination_number_offsets = [int(cn*(cn + 1)/2) for cn in range(13)]
        self.element_offset = 91
        
    def create(self, particle):
        system = Nanoparticle.from_atoms(particle)
        connectivity_matrix = system.get_connectivity_matrix()
        oc_a, oc_b = system.get_symbols_lists()
        ext_etop = EnvironmentalTopologiesExt(connectivity_matrix, oc_a, oc_b)
        ext_etop.compute_local_environments()

specs = [
    ('c_matrix', numba.int16[:,:]),
    ('oc_a', numba.float64[:]),
    ('oc_b', numba.float64[:])
]

#@jitclass(specs)
class EnvironmentalTopologiesExt():
    def __init__(self, c_matrix, oc_a, oc_b):
        self.c_matrix = c_matrix
        self.oc_a = oc_a
        self.oc_b = oc_b
        

    def compute_local_environments(self):

        for i in range(self.c_matrix.shape[0]):
            a_bonds = np.sum(self.c_matrix[i]*self.oc_a)
            b_bonds = np.sum(self.c_matrix[i]*self.oc_a)

            #self.atom_features[i][0] = a_bonds
            #self.atom_features[i][1] = b_bonds

    
    def compute_local_environment(self, system, atom_index):

        a_bonds = np.int16(np.sum(system.connectivity_matrix[atom_index]*system.occupancy_symbol_a))
        b_bonds = np.int16(np.sum(system.connectivity_matrix[atom_index]*system.occupancy_symbol_b))

        self.local_environments[atom_index][0] = a_bonds
        self.local_environments[atom_index][1] = b_bonds
        
    
    def compute_atom_features(self, system):

        for atom_idx in range(self.n_atoms):

            element = np.int16(system.occupancy_symbol_a[atom_idx])
            cn = np.sum(self.local_environments[atom_idx],dtype=int)
            off_set = self.element_offset * element
            self.atom_features[atom_idx] = np.sum([off_set + cn + self.local_environments[atom_idx][0]])
