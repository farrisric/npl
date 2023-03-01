from typing import List

import numpy as np

from sklearn import linear_model

from npl.calculators.base_calculator import BaseCalculator
from npl.core.nanoparticle import Nanoparticle


class TrainedCalculator(BaseCalculator):
    """A LinearModelCalculator object that calculates the 
    total energy of a nanoparticle based on its descriptor.

    Parameters
    ----------
    training_set : npl.Nanoparticle
        training set data of Nanoparticle objects that contains 
        the descriptors and the energ specified in energy_key.
    descriptor_key : str
        identifier to the descriptor used to train the linear model
    linear_model : sklearn.linear_model.*
        linear model used to fit the training_set data
    name : str
        identifier of the calculator

    """

    def __init__(self, name: str, linear_model) -> None:
        super().__init__(name=name)
        self.model = linear_model()
       
    
    def fit(self, training_set : List[Nanoparticle]) -> None:
        X, y = [], []
        for particle in training_set:
            X.append(particle.descriptors[self.name])
            y.append(particle.get_potential_energy())
        self.model.fit(X,y)

    def get_coefficients(self) -> np.array:
        return self.model.coef_

    def calculate_total(self, particle):
        return np.dot(np.transpose(self.model.coef_),particle.descriptors[self.name])
