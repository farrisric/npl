import numpy as np
from collections import defaultdict

from ase import Atoms
from ase.data import atomic_numbers, chemical_symbols

# TODO Derive directly from Atoms


class AtomWrapper:
    """Wrapper class around the ase Atoms object to extend it with some higher-level functionality.

    This class extends an Atoms object with additional functions for ordering, stoichiometry
    handling and is split from BaseNanoparticle to keep the size of the modules as small as
    possible.
    """
    def __init__(self):
        self.atoms = Atoms()

        return

    def add_atoms(self, new_atoms):
        self.atoms.extend(new_atoms)

    def get_atoms(self, indices):
        return self.atoms[indices]

    def get_position(self, idx):
        return self.atoms.positions[idx]

    def remove_atoms(self, indices):
        del self.atoms[indices]

    def get_ase_atoms(self, indices=None):
        """Return atoms specified by the indices.

        Parameters:
            indices: list/array of int
        """
        if indices is None:
            return self.atoms
        return self.atoms[indices]

    def swap_symbols(self, index_pairs):
        numbers = self.atoms.numbers
        for idx1, idx2 in index_pairs:
            numbers[idx1], numbers[idx2] = numbers[idx2], numbers[idx1]

    def random_ordering(self, new_stoichiometry):
        """Creates a random chemical ordering for the given stoichiometry

        The new stoichiometry has to be given as dict of type str -> int that matches
        the symbol to the number of atoms.

        Parameters:
            new_stoichiometry: dict
        """
        new_symbols = []
        for symbol in new_stoichiometry:
            new_symbols += [symbol]*new_stoichiometry[symbol]
        np.random.shuffle(new_symbols)
        self.atoms.symbols = new_symbols

    def transform_atoms(self, atom_indices, new_symbols):
        numbers = self.atoms.numbers
        for idx, symbol in zip(atom_indices, new_symbols):
            numbers[idx] = atomic_numbers[symbol]

    def get_indices(self):
        """Convenience function for range(n_atoms)."""
        return np.arange(0, len(self.atoms))

    def get_all_symbols(self):
        """Return list of symbols which occur at least once in the particle."""
        return list(self.atoms.symbols.species())

    def get_symbol(self, atom_idx):
        return chemical_symbols[self.atoms.numbers[atom_idx]]

    def get_symbols(self, indices=None):
        if indices is None:
            return self.atoms.symbols
        return self.atoms[indices].symbols

    def get_indices_by_symbol(self, symbol):
        """Return atom indices of given symbol.

        Parameters:
            symbol : str
        """
        indices = self.atoms.symbols.indices()
        return indices.get(symbol, np.array([]))

    def get_indices_by_symbol_map(self):
        """Return ``{symbol: indices}`` for every present symbol in a single pass.

        Cheaper than calling ``get_indices_by_symbol`` per symbol, which rebuilds
        the whole map each time.
        """
        return self.atoms.symbols.indices()

    def get_n_atoms(self):
        """Return the number of atoms."""
        return len(self.atoms)

    # TODO implement more efficiently using in-built ase methods
    def get_n_atoms_of_symbol(self, symbol):
        return len(self.get_indices_by_symbol(symbol))

    def get_stoichiometry(self):
        """Return the current composition as defaultdict.

        For symbols that are not in the particle, the defaultdict will evaluate to 0.
        """
        stoichiometry = self.atoms.symbols.indices()
        for symbol in stoichiometry:
            stoichiometry[symbol] = len(stoichiometry[symbol])
        return defaultdict(lambda: 0, stoichiometry)

    def get_positions(self, indices=None):
        """Return the positions of the atoms.

        If indices is not specified, all atoms will be returned.

        Parameters:
            indices : list/array of int
        """
        if indices is None:
            return self.atoms.positions
        return self.atoms[indices].positions

    def translate_atoms_positions(self, position):
        self.atoms.translate(-position)
