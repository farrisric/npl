from abc import ABC, abstractmethod

from npl.descriptors.base_descriptor import BaseDescriptor

class BaseSurrogateEnergyModel(ABC):
    """An abstract base class for creating surrogate energy models.
    
    Parameters
    ----------
    name:
        identifier of the descriptor

    """

    def __init__(self, name = 'BaseSurrogateEnergyModel'):
        self.name = name
        self._descriptor = BaseDescriptor
        
    @abstractmethod
    def train(self):
        pass

    @property
    def descriptor(self):
        return self._descriptor
    
        