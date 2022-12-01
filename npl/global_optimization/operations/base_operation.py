from abc import ABC, abstractmethod

class BaseOperator(ABC):
    """An abstract base class for operation that can be performed 
    in a global optimization run.
    
    Parameters
    ----------
    name:
        identifier of the descriptor
    """

    def __init__(self):
        pass

    @abstractmethod
    def perform_operation(self):
        pass

    @abstractmethod
    def revert_operation(self):
        pass