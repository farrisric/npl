from ase.neighborlist import natural_cutoffs
from ase.neighborlist import build_neighbor_list

from collections import defaultdict
import numpy as np
import itertools
import copy


def construct_neighborlist( atoms, scale_factor=1.0, npl = True, cutoffs = None):

        neighborlist = defaultdict(lambda: set())

        if not cutoffs:
            cutoffs = natural_cutoffs(atoms, mult=scale_factor)
            
        neighbor_list = build_neighbor_list(atoms,
                                            cutoffs=cutoffs,
                                            bothways=True,
                                            self_interaction=False)

        for atom_idx, _ in enumerate(atoms):
            neighbors, _ = neighbor_list.get_neighbors(atom_idx)
            if npl:
                neighborlist[atom_idx] = set(neighbors)
            else:
                neighborlist[atom_idx] = neighbors
        return neighborlist
