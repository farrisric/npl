from sklearn import linear_model

from npl.calculators.base_calculator import BaseCalculator
from npl.core.nanoparticle import Nanoparticle

class LinearModelCalculator(BaseCalculator):
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

    def __init__(self,
                training_set: list[Nanoparticle],
                descriptor_key: str,
                linear_model: linear_model,
                name: str):

        super().__init__(name=name)

        self.model = linear_model
        self.descriptor_key = descriptor_key
        self._fit(training_set, descriptor_key)        
    
    def _fit(self, trainin_set, descriptor_key):

        X, y = [], []
        for particle in trainin_set:
            X.append(particle.info['descriptors'][descriptor_key])
            y.append(particle.get_potential_energy())

        self.model.fit(X,y)

    def calculate_total(self, particle):
        return self.model.predict(particle.info['descriptors'][self.descriptor_key])
