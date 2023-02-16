from npl.core import Nanoparticle
from npl.global_optimization.operations import BaseOperator
from sortedcontainers import SortedKeyList
from typing import List
import numpy as np

class GuidedExchangeOperator(BaseOperator):
    def __init__(self, system: Nanoparticle, environment_energies : List[float]) -> None:
        super().__init__(system)
        self.n_envs = int(len(environment_energies)/len(self.atomic_numbers))
        self.flip_energies = {atomic_number : None for atomic_number in self.atomic_numbers}
        
        self.get_flip_energies(environment_energies)

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
        
    
    
