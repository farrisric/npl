import numpy as np
from numba import njit

from npl.core import Nanoparticle
from npl.descriptors import BaseDescriptor


class EnvironmentalTopologies(BaseDescriptor):
    def __init__(self):
        super().__init__(name='ETOP')
        

    def create(self, particle):
        system = Nanoparticle.from_atoms(particle)
        bond_matrix = system.get_bond_matrix()
        feature_vector = self.compute_feature_vector(system._connectivity_matrix, system._occupation_matrix.T, bond_matrix)
        return feature_vector

    def _create(self, particle):
        system = Nanoparticle.from_atoms(particle)
        bond_matrix = system.get_bond_matrix()
        feature_vector = self.compute_feature_vector(system._connectivity_matrix, system._occupation_matrix.T, bond_matrix)
        particle.info['descriptors'] = {self.name : feature_vector}

    @staticmethod
    @njit
    def compute_feature_vector(occupation_matrix, connectivity_matrix, bond_matrix):
        coordination_number_offsets = np.array([int(cn*(cn + 1)/2) for cn in range(13)])
        element_offset = 91
        feature_vector = np.zeros(182)  
        for symbol, bonds in zip(occupation_matrix[1], bond_matrix):        
            index = int(coordination_number_offsets[int(sum(bonds))] + bonds[0] + (element_offset*symbol))
            feature_vector[index] += 1
        return feature_vector

    def update(self, system: Nanoparticle):
        bond_matrix = system.get_bond_matrix()
        feature_vector = self.compute_feature_vector(system._connectivity_matrix, system._occupation_matrix.T, bond_matrix)
        return feature_vector

    
