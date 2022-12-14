from abc import ABC, abstractmethod


class BaseCalculator(ABC):
    """
    Base class for calculators.

    Attributes
    ----------
    name : str
        human-readable calculator name
    """

    def __init__(self, name='BaseCalculator'):
        self.name = name

    @abstractmethod
    def calculate_total(self):
        pass

    @abstractmethod
    def calculate_change(self):
        pass

    @property
    def sublattices(self):
        raise NotImplementedError()
