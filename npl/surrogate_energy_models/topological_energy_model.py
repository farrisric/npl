from typing import List

#from sklearn.linear_model
import numpy as np
from ase import Atoms

from npl.core import Nanoparticle 
from npl.surrogate_energy_models import BaseSurrogateEnergyModel
from npl.descriptors import Topologies

class TopologicalEnergyModel(BaseSurrogateEnergyModel):
    """ A class to create a Topological EnergyModel
    """

    def __init__(self, symbols : List[str]) -> None:
        super().__init__(name = 'Topological Energy Model')
        self._descriptor = Topologies(symbols)

    def get_training_set(self, training_set : List[Atoms]) -> List[np.array]:

        X = []
        for atoms in training_set:
            nanoparticle = Nanoparticle.from_atoms(atoms)
            X.append(self.descriptor.create(nanoparticle))

        return X

    def train(self, training_set : List[Atoms]):
        return super().train()