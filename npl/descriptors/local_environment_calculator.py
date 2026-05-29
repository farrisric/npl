import numpy as np
import copy


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

    def predict_local_environment(self, particle, lattice_index):
        n_a_atoms = 0
        n_b_atoms = 0

        neighbors = particle.neighbor_list[lattice_index]
        for neighbor in neighbors:
            if particle.get_symbol(neighbor) == self.symbol_a:
                n_a_atoms += 1
            else:
                n_b_atoms += 1

        return np.array([n_a_atoms, n_b_atoms])
