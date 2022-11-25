from numba import njit
import numpy as np

from npl.descriptors import Descriptor
from npl.core import Nanoparticle

class Topologies(Descriptor):
    def __init__(self):
        pass

    def create(self, particle):
        system = Nanoparticle.from_atoms(particle)
        bond_matrix = system.get_bond_matrix()
        feature_vector = self.compute_feature_vector(system._connectivity_matrix, system._occupation_matrix.T, bond_matrix)
        return feature_vector

    def __init__(self):
        self.neighbor_list = None
        self.symbols = None
        self.symbols_list = None
        super().__init__()

    def create(self, particle):

        self.get_connectivity_matrix(particle)
        self.get_symbols_list(particle)
        #return super().create(particle)

    def get_symbols_list(self, particle):
        self.symbols = sorted(list(particle.symbols.indices()))
        self.symbols_list = particle.symbols

    def get_connectivity_matrix(self,particle):
        self.neighbor_list = get_connectivity_matrix(particle)
    
    def count_bonds(self, particle):
        n_aa_bonds = 0
        n_ab_bonds = 0
        n_bb_bonds = 0

        for lattice_index_with_symbol_a in particle.atoms.get_indices_by_symbol(self.symbols[0]):
            neighbor_list = self.neighbor_list[lattice_index_with_symbol_a]
            for neighbor in neighbor_list:
                symbol_neighbor = particle.atoms.get_symbol(neighbor)

                if self.symbol_a != symbol_neighbor:
                    n_ab_bonds += 0.5
                else:
                    n_aa_bonds += 0.5

        for lattice_index_with_symbol_b in particle.atoms.get_indices_by_symbol(self.symbols[1]):
            neighbor_list = self.neighbor_list[lattice_index_with_symbol_b]
            for neighbor in neighbor_list:
                symbol_neighbor = particle.atoms.get_symbol(neighbor)

                if self.symbol_b == symbol_neighbor:
                    n_bb_bonds += 0.5
                else:
                    n_ab_bonds += 0.5

        return n_aa_bonds, n_bb_bonds, n_ab_bonds

    def count_coordination_occurrences():
        pass

