from ase.neighborlist import natural_cutoffs
from ase.neighborlist import build_neighbor_list

from collections import defaultdict
import numpy as np
import itertools
import copy

def get_connectivity_matrix(atoms, scale_factor=1.0, cutoffs = None):

        connectivity_matrix = np.zeros((len(atoms), len(atoms)), dtype=np.int128)
            
        neighbor_list = build_neighbor_list(atoms,
                                            cutoffs=cutoffs,
                                            bothways=True,
                                            self_interaction=False)

        for atom_idx, _ in enumerate(atoms):
            neighbors, _ = neighbor_list.get_neighbors(atom_idx)
            for neighbor_idx in neighbors:
                connectivity_matrix[atom_idx][neighbor_idx] = 1

        return connectivity_matrix


