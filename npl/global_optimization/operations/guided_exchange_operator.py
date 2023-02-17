from npl.core import Nanoparticle
from npl.global_optimization.operations import BaseOperator
from npl.global_optimization.operations import features
from sortedcontainers import SortedKeyList
from typing import List, Dict
import numpy as np

class GuidedExchangeOperator(BaseOperator):
    def __init__(self, system: Nanoparticle, environment_energies : List[float], feature_index_values : Dict[tuple, float]) -> None:
        super().__init__(system)
        self.feature_index_values = feature_index_values
        self.exchange_energies = dict()
        self.sorted_indices = SortedKeyList(key=lambda x: min(self.exchange_energies[x]))
        self.n_envs = int(len(environment_energies)/len(self.atomic_numbers))
        self.flip_energies = {atomic_number : None for atomic_number in self.atomic_numbers}

        self.get_flip_energies(environment_energies)
        self.bind_system(system)

    def bind_system(self, system):
        for atom in system:
            Z = atom.number
            number_index = system.get_number_index(Z)
            atom_feature = self.get_atom_feature(system, atom)
            feature_index = self.feature_index_values[Z][atom_feature]
            self.exchange_energies[atom.index] = self.flip_energies[Z][feature_index]
            self.sorted_indices.add(atom.index)

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

    
    def revert_operation(self):
        return super().revert_operation()

    def perform_operation(self):
        return super().perform_operation()
        
    
    
