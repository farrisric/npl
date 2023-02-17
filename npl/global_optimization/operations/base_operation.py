from abc import ABC, abstractmethod
from itertools import combinations
from typing import List

from npl.core import Nanoparticle


class BaseOperator(ABC):
    """An abstract base class for operation that can be performed 
    in a global optimization run.
    
    Parameters
    ----------
    name:
        identifier of the descriptor
    """

    def __init__(self, system : Nanoparticle) -> None:
        self.exchange_types = []
        self.indices_by_element = dict()
        self.atomic_numbers = system.get_unique_atomic_numbers()

        self.get_indices_by_symbol(system)
        self.compute_number_of_exchange_types(system)

    @abstractmethod
    def perform_operation(self):
        pass

    @abstractmethod
    def revert_operation(self):
        pass
    
    # def bind_system(self, system):
    #     """Bind the nanoparticle to the exchange operator in 
    #     order to perform swaps
    #     """
    #     self.get_indices_by_symbol(system)
    #     self.compute_number_of_exchange_types(system)


    def compute_number_of_exchange_types(self, system):
        """Compute number of possible exchanges it can be performs.
        Depends on the number of unqiue elements in the system:
        e.g. 3 elements generate 3 type of possible exchanges
        """
        self.exchange_types =  [ex_type for ex_type in combinations(range(len(self.indices_by_element)) ,2)]

    def get_indices_by_symbol(self, system):
        """Separates in lists the atom indices for each element in the system
        """
        for atomic_number in self.atomic_numbers:
            element_indices = system.get_occupation_indices_by_symbol(atomic_number)
            number_index = system.get_number_index(atomic_number)
            self.indices_by_element[number_index] = element_indices
    
    def execute_swap_operation(self, a, exchange, swap):
        ex1, ex2 = exchange
        swap1, swap2 = swap
        a[ex1][swap1], a[ex2][swap2] = a[ex2][swap2], a[ex1][swap1]