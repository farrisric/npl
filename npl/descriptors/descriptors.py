from abc import ABC, abstractmethod

class Descriptor(ABC):
    """An abstract base class for all descriptors."""

    def __init__(self):
        pass

    @abstractmethod
    def create(self, particle):
        """Creates the descriptor for the given systems.

        Args:
            particle (ase.Atoms): The nanoparticle for which to create the descriptor.
        """