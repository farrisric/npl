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
        self.sorted_indices = SortedKeyList(key=lambda x: min(self.exchange_energies[x]))
        self.n_envs = int(len(environment_energies)/len(self.atomic_numbers))
        self.flip_energies = {atomic_number : None for atomic_number in self.atomic_numbers}
        self.sorted_energies_dict = {Z : {} for Z in range(self.n_symbols)}
        self.swap_combinations = [x for x in get_combinations(range(self.n_symbols), 2) if x[1]!=x[0]] 
        for Z_i in self.sorted_energies_dict:
            for Z_j in self.sorted_energies_dict:
                if Z_j == Z_i:
                    continue
                self.sorted_energies_dict[Z_i][Z_j] = SortedKeyList(key=lambda x: self.exchange_energies[x][Z_j])

        self.get_flip_energies(environment_energies)
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
        atom_feature = self.get_atom_feature(system, atom)
        feature_index = self.feature_index_values[atom.number][atom_feature]
        self.exchange_energies[atom.index] = self.flip_energies[atom.number][feature_index]

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
            else: 
                return False
        return exchange_indices

    def get_atom_feature(self, system, atom):
        from npl.global_optimization.operations import features
        n = len(self.atomic_numbers)
        i = features[n][system.get_number_index(atom.number)]
        return tuple(system.atom_features['TOP'][atom.index][i])
         
    def env_from_symbol(self, index : int) -> float:
        return self.n_envs * index

    def get_flip_energies(self, environment_energies : List[float]):
        for Z in self.atomic_numbers:
            flip_energy_matrix = np.empty((self.n_envs, len(self.atomic_numbers)))
            for i in range(self.n_envs):
                for j in range(len(self.atomic_numbers)):
                    flip_initial = i + self.env_from_symbol(self.atomic_numbers.index(Z))
                    flip_final = i + self.env_from_symbol(j)
                    flip = environment_energies[flip_final] - environment_energies[flip_initial]
                    flip_energy_matrix[i][j] = flip
            self.flip_energies[Z] = flip_energy_matrix

#TODO impolement basin hopping steps and update function after swap

    def basin_hop_step(self):
        raise NotImplementedError
    
    def update(self, system : Nanoparticle, exchaged_indices : List[int], old_numbers : List[int]) -> None:
        for index, old_num in zip(exchaged_indices, old_numbers):
            Z_i = system.get_index_number(old_num)
            for Z_j in range(self.n_symbols):
                if Z_j != Z_i:
                    self.sorted_energies_dict[Z_i][Z_j].remove(index)
            self.get_exchange_energy(system, system[index])
            self.add_atom_index(system, system[index])
        
    def revert_operation(self):
        return super().revert_operation()

    def perform_operation(self, system : Nanoparticle):
        exchaged_indices = self.find_best_swap_pair()
        if not exchaged_indices:
            return print('Cannot do more steps')
        old_numbers = [system[x].number for x in exchaged_indices]
        i, j = exchaged_indices
        system[i].number, system[j].number = system[j].number, system[i].number
        system.update_occupation_matrix(exchaged_indices, old_numbers)
        self.update(system, exchaged_indices, old_numbers)

        
    
    
