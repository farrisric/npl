import numpy as np
import copy

from ase.data import atomic_numbers


class LocalEnvironmentCalculator:
    """Base class for local environment calculators.

    Intended use: Implementations should implement predict_local_environment() which calculates
    and returns the local environment of a single atom.
    """
    def __init__(self):
        pass

    def compute_local_environments(self, particle):
        for index in particle.get_indices():
            self.compute_local_environment(particle, index)

    def compute_local_environment(self, particle, atom_index):
        local_env = self.predict_local_environment(particle, atom_index)
        particle.set_local_environment(atom_index, local_env)

    def predict_local_environment(self, particle, lattice_index):
        raise NotImplementedError


class NeighborCountingEnvironmentCalculator(LocalEnvironmentCalculator):
    """Calculate a local environment of the form [n_a, n_b], where n_a denotes the number of
    surrounding a atoms.

    Currently restricted to two elements which are ordered alphabetically. Requires a valid
    neighbor list.
    """
    def __init__(self, symbols):
        LocalEnvironmentCalculator.__init__(self)
        symbols_copy = copy.deepcopy(symbols)
        symbols_copy.sort()
        self.symbol_a = symbols_copy[0]
        self.symbol_b = symbols_copy[1]
        self.number_a = atomic_numbers[self.symbol_a]

    def predict_local_environment(self, particle, lattice_index):
        # Count symbol_a neighbors directly on the atomic-numbers array; the rest
        # are symbol_b (this calculator is restricted to two elements).
        numbers = particle.atoms.atoms.numbers
        neighbors = particle.neighbor_list[lattice_index]
        n_a_atoms = sum(1 for neighbor in neighbors if numbers[neighbor] == self.number_a)
        return np.array([n_a_atoms, len(neighbors) - n_a_atoms])
