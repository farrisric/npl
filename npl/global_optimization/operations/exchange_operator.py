from itertools import combinations
import numpy as np

from npl.global_optimization.operations.base_operation import BaseOperator

class ExchangeOperator(BaseOperator):
    """Class that performs exchange operators between atoms of
    different elements in a nanoparticle. It uses the Occupation Matrix
    of a system to perform element swaps
    
    Parameters
    ----------
    name:
        identifier of the descriptor
    """

    def __init__(self, system, p_geometric):
        self.exchange_types = []
        self.indices_by_element = dict()
        self.p_geometric = p_geometric
        self.operations = []
        
        self.bind_system(system)

        
    def _compute_number_of_exchange_types(self, system):
        """Compute number of possible exchanges it can be performs.
        Depends on the number of unqiue elements in the system:
        e.g. 3 elements generate 3 type of possible exchanges
        """
        self.exchange_types =  [ex_type for ex_type in combinations(range(len(self.indices_by_element)) ,2)]

    def _get_indices_by_symbol(self, system):
        """Separates in lists the atom indices for each element in the system
        """
        for atomic_number in system.get_unique_numbers():
            element_indices = system.get_occupation_indices_by_symbol(atomic_number)
            number_index = system.get_number_index(atomic_number)
            self.indices_by_element[number_index] = element_indices

    def bind_system(self, system):
        """Bind the nanoparticle to the exchange operator in 
        order to perform swaps
        """
        self._get_indices_by_symbol(system)
        self._compute_number_of_exchange_types(system)

    def get_random(self, a, size = 1, replace=False):
        return np.random.choice(a, size, replace=replace)

    def get_random_exchange(self):
        n_exchanges = np.random.geometric(p=self.p_geometric, size=1)[0]
        exchanges = self.get_random(len(self.exchange_types), size = n_exchanges, replace=True)
        exchanges = [self.exchange_types[x] for x in exchanges]
        return exchanges

    def get_random_swap(self, exchage):

        symbol_index1, symbol_index2 = exchage
        idx1 = self.get_random(len(self.indices_by_element[symbol_index1]))[0]
        idx2 = self.get_random(len(self.indices_by_element[symbol_index2]))[0]

        swap_OM = (self.indices_by_element[symbol_index1][idx1], self.indices_by_element[symbol_index2][idx2])
        swap_indices_list = (idx1, idx2)
        return swap_OM, swap_indices_list

    def execute_swap_operation(self, a, exchange, swap):
        ex1, ex2 = exchange
        swap1, swap2 = swap
        
        a[ex1][swap1], a[ex2][swap2] = a[ex2][swap2], a[ex1][swap1]


    def perform_operation(self, system):
        """Perform the swap operation for a Monte Carlo step"""
        exchanges = self.get_random_exchange()
        for exchange in exchanges:
            swap_OM, swap_indices_list = self.get_random_swap(exchange)
            self.operations.append((exchange, swap_OM, swap_indices_list))

            self.execute_swap_operation(system.OM, exchange, (swap_OM[0], swap_OM[0]))
            self.execute_swap_operation(system.OM, exchange, (swap_OM[1], swap_OM[1]))
            self.execute_swap_operation(self.indices_by_element, exchange, swap_indices_list)

    def revert_operation(self, system):
        """Revert the swap operation of a Monte Carlo step"""
        for operation in self.operations:
            exchange, swap_OM, swap_indices_list = operation

            self.execute_swap_operation(system.OM, exchange, (swap_OM[0], swap_OM[0]))
            self.execute_swap_operation(system.OM, exchange, (swap_OM[1], swap_OM[1]))
            self.execute_swap_operation(self.indices_by_element, exchange, swap_indices_list)


    # def get_exchange_types(self):
    #     n_exchanges = np.random.geometric(p=self.p_geometric, size=1)[0]
    #     exchange_idx = np.random.choice(len(self.exchange_types), n_exchanges)
    #     exchange_types = [self.exchange_types[x] for x in exchange_idx]
    #     return exchange_types

    # def get_swap_indices(self):
    #     exchange_types = self.get_exchange_types()
    #     swaps = []
    #     for exchange_type in exchange_types:
    #         symbol1, symbol2 = exchange_type
    #         symbol1_indices = np.random.choice(self.indices_by_element[symbol1], replace=False)
    #         symbol2_indices = np.random.choice(self.indices_by_element[symbol2], replace=False)

    #         swaps.append((symbol1_indices,symbol2_indices))
    #     return swaps, exchange_types

    # def get_operation(self):
    #     swaps, exchange_types = self.get_swap_indices()
    #     self.operation = swaps, exchange_types
    #     return swaps, exchange_types

    # def exchange_symbols(self, system, swaps, exchange_types):
    #     for swap, exchange_type in zip(swaps, exchange_types):
    #         idx1, idx2 = swap
    #         symbol1, symbol2 = exchange_type
    #         system.OM[symbol1][idx1], system.OM[symbol2][idx2] = system.OM[symbol1][idx2], system.OM[symbol2][idx1]

    # def _update_indices(self, swaps, exchange_types):

    #     self.indices_by_element


    # def perform_operation(self, system):
    #     swaps, exchange_types = self.get_operation()
    #     self.exchange_symbols(system, swaps, exchange_types)

    # def revert_operation(self, system):
    #     swaps, exchange_types = self.operation
    #     self.exchange_symbols(system, swaps, exchange_types)

