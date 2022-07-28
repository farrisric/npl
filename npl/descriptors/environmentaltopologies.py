from numba import njit
import numpy as np

from descriptors import Descriptor
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

    @njit
    def compute_local_environments(self, system):

        for i in range(len(system.connectivity_matrix)):
            a_bonds = np.int16(np.sum(system.connectivity_matrix[i]*system.occupancy_symbol_a))
            b_bonds = np.int16(np.sum(system.connectivity_matrix[i]*system.occupancy_symbol_b))

            self.atom_features[i][0] = a_bonds
            self.atom_features[i][1] = b_bonds

    @njit
    def compute_local_environment(system, atom_index):

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
