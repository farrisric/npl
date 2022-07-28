from ase.neighborlist import natural_cutoffs
from ase.neighborlist import build_neighbor_list

from collections import defaultdict
import numpy as np
import itertools
import copy


def get_connectivity_matrix(system, scale_factor=1.0, cutoffs = None):

        connectivity_matrix = np.zeros((len(system), len(system)), dtype=np.int16)
            
        neighbor_list = build_neighbor_list(system,
                                            cutoffs=cutoffs,
                                            bothways=True,
                                            self_interaction=False)

        for atom_idx, _ in enumerate(system):
            neighbors, _ = neighbor_list.get_neighbors(atom_idx)
            for neighbor_idx in neighbors:
                connectivity_matrix[atom_idx][neighbor_idx] = 1

        return connectivity_matrix


def get_symbols_lists(system):
    position_symbol_a = np.zeros(len(system))
    position_symbol_b = np.zeros(len(system))

    for atom in system:
        if atom.symbol == list(system.symbols.indices())[0]:
            position_symbol_a[atom.index] = 1
        else:
            position_symbol_b[atom.index] = 1
