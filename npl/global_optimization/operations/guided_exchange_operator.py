from sortedcontainers import SortedKeyList
import numpy as np

class GuidedExchangeOperator:
    def __init__(self, environment_energies, symbols):
        self.symbols = sorted(symbols)
        self.n_envs = int(len(environment_energies)/len(symbols))


    def get_flip_energies(self, environment_energies):

        flip_energies_by_symbols = {symbol : np.empty((len(self.symbols), self.n_envs)) for symbol in self.symbols}
        
        a_off = 0
        for symbol in self.symbols:
            
            for i in range(self.n_envs):
                for j in range(len(self.symbols)):

                    flip_energies_by_symbols[symbol][i][j] = environment_energies[i] - environment_energies[i]
        



        self.env_energy_differences = [environment_energies[i] - environment_energies[i + self.n_envs] for i in range(self.n_envs)]


        self.symbol1_exchange_energies = dict()
        self.symbol2_exchange_energies = dict()

        self.symbol1_indices = SortedKeyList(key=lambda x: self.symbol1_exchange_energies[x])
        self.symbol2_indices = SortedKeyList(key=lambda x: self.symbol2_exchange_energies[x])

        self.n_symbol1_atoms = 0
        self.n_symbol2_atoms = 0