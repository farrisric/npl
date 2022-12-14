from typing import List

from numba import njit
import numpy as np
from itertools import combinations_with_replacement

from ase.data import atomic_numbers

from npl.descriptors import BaseDescriptor
from npl.core import Nanoparticle

class Topologies(BaseDescriptor):
    def __init__(self, symbols : List[str]) -> None:
        self.n_features = 0

        self.atomic_numbers = sorted([atomic_numbers[symbol] for symbol in symbols])
        self.bond_types = None
        self.bond_type_per_element = None

        super().__init__(name = 'TOP')
        self._compute_bond_types()
        self._compute_n_features()

    def _compute_n_features(self):
        self.n_features += len(self.bond_types)
        self.n_features += len(self.atomic_numbers) * 13

    def _compute_bond_types(self):

        self.bond_types = [x for x in combinations_with_replacement(self.atomic_numbers, 2)]
        self.bond_type_per_element = {x : [] for x in self.atomic_numbers}

        for atomic_number in self.atomic_numbers:
            for bond_type in self.bond_types:
                if atomic_number in bond_type:
                    self.bond_type_per_element[atomic_number].append(self.bond_types.index(bond_type))
    
    def _compute_atom_features(self, system, atom_index, row):

        key = system[atom_index].number

        for element, feature_index in enumerate(self.bond_type_per_element[key]):
            row[feature_index] += np.dot(system.CM[atom_index], system.OM[element])

    def _compute_all_atom_features(self, system):
        system.atom_features[self.name] = np.zeros((len(system), self.n_features), dtype=np.int8)
  
        for atom_index in range(len(system)):
            self._compute_atom_features(system, atom_index, system.atom_features[self.name][atom_index])

    def get_feature_vector(self, system : Nanoparticle):

        if not isinstance(system, Nanoparticle):
            system = Nanoparticle.from_atoms(system)

        self._compute_all_atom_features(system)
        return system.atom_features[self.name].sum(axis=0)
        
        


    


    # def create(self, particle):
    #     system = Nanoparticle.from_atoms(particle)
    #     bond_matrix = system.get_bond_matrix()
    #     feature_vector = self.compute_feature_vector(system._connectivity_matrix, system._occupation_matrix.T, bond_matrix)
    #     return feature_vector

    

    # def create(self, particle):

    #     self.get_connectivity_matrix(particle)
    #     self.get_symbols_list(particle)
    #     #return super().create(particle)

    # def get_symbols_list(self, particle):
    #     self.symbols = sorted(list(particle.symbols.indices()))
    #     self.symbols_list = particle.symbols

    # def get_connectivity_matrix(self,particle):
    #     self.neighbor_list = get_connectivity_matrix(particle)
    
    # def count_bonds(self, particle):
    #     n_aa_bonds = 0
    #     n_ab_bonds = 0
    #     n_bb_bonds = 0

    #     for lattice_index_with_symbol_a in particle.atoms.get_indices_by_symbol(self.symbols[0]):
    #         neighbor_list = self.neighbor_list[lattice_index_with_symbol_a]
    #         for neighbor in neighbor_list:
    #             symbol_neighbor = particle.atoms.get_symbol(neighbor)

    #             if self.symbol_a != symbol_neighbor:
    #                 n_ab_bonds += 0.5
    #             else:
    #                 n_aa_bonds += 0.5

    #     for lattice_index_with_symbol_b in particle.atoms.get_indices_by_symbol(self.symbols[1]):
    #         neighbor_list = self.neighbor_list[lattice_index_with_symbol_b]
    #         for neighbor in neighbor_list:
    #             symbol_neighbor = particle.atoms.get_symbol(neighbor)

    #             if self.symbol_b == symbol_neighbor:
    #                 n_bb_bonds += 0.5
    #             else:
    #                 n_ab_bonds += 0.5

    #     return n_aa_bonds, n_bb_bonds, n_ab_bonds

    # def count_coordination_occurrences():
    #     pass

