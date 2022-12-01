from abc import ABC, abstractmethod

class Descriptor(ABC):
    """An abstract base class for all descriptors.
    
    Parameters
    ----------
    name:
        identifier of the descriptor
    """

    def __init__(self, name):
        self.name = name
        self.n_features = 0

    @abstractmethod
    def _compute_n_features(self):
        pass

    @abstractmethod
    def create(self, particle):
        """Creates the descriptor for the given systems.

        Args:
            particle (ase.Atoms): The nanoparticle for which to create the descriptor.
        """

    