from typing import Dict, List

import numpy as np
from sortedcontainers import SortedKeyList

from npl.core import Nanoparticle
from npl.global_optimization.operations import BaseOperator, features
from npl.utils import get_combinations


class GuidedExchangeOperator(BaseOperator):
    def __init__(self, system : Nanoparticle, environment_energies : List[float], feature_index_values : Dict[tuple, float]) -> None:
        super().__init__(system)
        self.feature_index_values = feature_index_values
        self.exchange_energies = dict()
        self.n_envs = int(len(environment_energies)/len(self.atomic_numbers))
        self.env_energy_difference = None
        self.swap_combinations = [x for x in get_combinations(range(self.n_symbols), 2) if x[1] != x[0]] 
        
        self.sorted_energies_dict = self.get_sorted_dict_by_n_symbols(self.n_symbols)
        self.get_env_energy_difference(environment_energies)
        self.bind_system(system)

        self.exchanged_indices = []

    #TODO change exchange energy with flip energy
    def bind_system(self, system):
        "Sorts each atom in the system based on each flip energy it has"
        for atom in system:
            self.get_exchange_energy(system, atom)
            self.add_atom_index(system, atom)

    def get_exchange_energy(self, system, atom):
        "Computes the flip energies of atom_i based on its Z and features (bonds)"
        feature_index = self.get_feature_index_from_atom_feature(system, atom)
        self.exchange_energies[atom.index] = self.env_energy_difference[feature_index]

    def get_feature_index_from_atom_feature(self, system, atom):
        atom_feature = self.get_atom_feature(system, atom)
        return self.feature_index_values[atom.number][atom_feature]

    def add_atom_index(self, system, atom):
        "Add the atom.index into each sorted list"
        Z_i = system.get_number_index(atom.number)
        for Z_j in range(self.n_symbols):
            if Z_j != Z_i:
                self.sorted_energies_dict[Z_i][Z_j].add(atom.index)

    def find_best_swap_pair(self):
        "Finds the best pair of atoms to swap that yields the maximum decrease in energy"
        best_swap = 0
        for Z_i, Z_j in self.swap_combinations:
            a = self.sorted_energies_dict[Z_i][Z_j][0]
            b = self.sorted_energies_dict[Z_j][Z_i][0]
            flip = self.exchange_energies[a][Z_j] + self.exchange_energies[b][Z_i]
            if flip < best_swap:
                exchange_indices = (a, b)
                best_swap = flip
        if best_swap < 0:
            return exchange_indices, best_swap
        return False, False

    def get_atom_feature(self, system, atom):
        from npl.global_optimization.operations import features
        n = len(self.atomic_numbers)
        i = features[n][system.get_number_index(atom.number)]
        return tuple(system.atom_features['TOP'][atom.index][i])
         
    def env_from_symbol(self, index : int) -> float:
        return self.n_envs * index

    def get_env_energy_difference(self, environment_energies : List[float]):
        self.env_energy_difference = np.zeros((len(environment_energies), len(self.atomic_numbers)))
        for Z in range(len(self.atomic_numbers)):
            for i in range(self.n_envs):
                for j in range(len(self.atomic_numbers)):
                    flip_initial = i + self.env_from_symbol(Z)
                    flip_final = i + self.env_from_symbol(j)
                    flip = environment_energies[flip_final] - environment_energies[flip_initial]
                    self.env_energy_difference[flip_initial][j] = flip

#TODO implement basin hopping steps and update function after swap

    def basin_hop_step(self):
        raise NotImplementedError
    
    def update(self, system : Nanoparticle, exchanged_indices : List[int], old_numbers : List[int]) -> None:
        for index, old_num in zip(exchanged_indices, old_numbers):
            Z_i = system.get_index_number(old_num)
            for Z_j in range(self.n_symbols):
                if Z_j != Z_i:
                    self.sorted_energies_dict[Z_i][Z_j].remove(index)
            self.get_exchange_energy(system, system[index])
            self.add_atom_index(system, system[index])
        
    def revert_operation(self):
        return super().revert_operation()

    def perform_operation(self, system : Nanoparticle):
        exchanged_indices, best_swap = self.find_best_swap_pair()
        if not exchanged_indices:
            return 'no more step bitch'
        neighbors = self.get_neighbor_list(system, exchanged_indices)
        neighbors = self.get_neighbor_list(system, neighbors)
        old_numbers = [system[x].number for x in neighbors]
        i, j = exchanged_indices
        system[i].symbol, system[j].symbol = system[j].symbol, system[i].symbol
        system.update_occupation_matrix(neighbors, old_numbers)
        self.update(system, neighbors, old_numbers)
        return best_swap

    def get_sorted_dict_by_n_symbols(self, n_symbols): 
        sorted_energies_dict_2 = {
            0: {1 : SortedKeyList(key=lambda i: self.exchange_energies[i][1])},
            1: {0 : SortedKeyList(key=lambda i: self.exchange_energies[i][0])}
        }
        sorted_energies_dict_3 = {
            0 : {
                1 : SortedKeyList(key=lambda i: self.exchange_energies[i][1]),
                2 : SortedKeyList(key=lambda i: self.exchange_energies[i][2])
            },
            1 : {
                0 : SortedKeyList(key=lambda i: self.exchange_energies[i][0]),
                2 : SortedKeyList(key=lambda i: self.exchange_energies[i][2])
            },
            2: {
                0 : SortedKeyList(key=lambda i: self.exchange_energies[i][0]),
                1 : SortedKeyList(key=lambda i: self.exchange_energies[i][1])
            }
        }
        sorted_energies_dict_4 = {
            0 : {
                1 : SortedKeyList(key=lambda i: self.exchange_energies[i][1]),
                2 : SortedKeyList(key=lambda i: self.exchange_energies[i][2]),
                3 : SortedKeyList(key=lambda i: self.exchange_energies[i][3])
            },
            1 : {
                0 : SortedKeyList(key=lambda i: self.exchange_energies[i][0]),
                2 : SortedKeyList(key=lambda i: self.exchange_energies[i][2]),
                3 : SortedKeyList(key=lambda i: self.exchange_energies[i][3])
            },
            2: {
                0 : SortedKeyList(key=lambda i: self.exchange_energies[i][0]),
                1 : SortedKeyList(key=lambda i: self.exchange_energies[i][1]),
                3 : SortedKeyList(key=lambda i: self.exchange_energies[i][3])
            },
            3 : {
                0 : SortedKeyList(key=lambda i: self.exchange_energies[i][0]),
                1 : SortedKeyList(key=lambda i: self.exchange_energies[i][1]),
                2 : SortedKeyList(key=lambda i: self.exchange_energies[i][2])
            }
        }

        if n_symbols == 2:
            return sorted_energies_dict_2
        if n_symbols == 3:
            return sorted_energies_dict_3
        if n_symbols == 4:
            return sorted_energies_dict_4


        
    
    
