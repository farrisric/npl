from itertools import combinations_with_replacement

from npl.global_optimization.operations.base_operation import BaseOperator

class ExchangeOperator:
    """Class that performs exchange operators between atoms of
    different elements in a nanoparticle. It uses the Occupation Matrix
    of a system to perform element swaps
    
    Parameters
    ----------
    name:
        identifier of the descriptor
    """

    def __init__(self, system):
        pass

    def bind_system(self, system):
        exchange_types = 
