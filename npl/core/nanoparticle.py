import sys
from ase import Atoms
import numpy as np

class Nanoparticle(Atoms):
    """Nanoparticle Object used within the package, inerithms from the
    ase.Atoms class"""

    def __init__(
        self,
        symbols=None,
        positions=None,
        numbers=None,
        tags=None,
        momenta=None,
        masses=None,
        magmoms=None,
        charges=None,
        scaled_positions=None,
        cell=None,
        pbc=None,
        celldisp=None,
        constraint=None,
        calculator=None,
        info=None
    ):

        super().__init__(
            symbols,
            positions,
            numbers,
            tags,
            momenta,
            masses,
            magmoms,
            charges,
            scaled_positions,
            cell,
            pbc,
            celldisp,
            constraint,
            calculator,
            info,
        )

        self._connectivity_matrix = None
        self.occupancy_symbol_a = None
        self.occupancy_symbol_b = None

    @staticmethod
    def from_atoms(atoms):
        """Creates a Nanoparticle object from ASE.Atoms object."""
        system = Nanoparticle(
            symbols=atoms.get_chemical_symbols(),
            positions=atoms.get_positions(),
            tags=atoms.get_tags(),
            momenta=atoms.get_momenta(),
            masses=atoms.get_masses(),
            magmoms=atoms.get_initial_magnetic_moments(),
            charges=atoms.get_initial_charges(),
            cell=atoms.get_cell(),
            pbc=atoms.get_pbc(),
            celldisp=atoms.get_celldisp(),
            constraint=atoms._get_constraints(),
            calculator=atoms.get_calculator(),
            info=atoms.info,
        )

        return system
    
    def get_symbols_lists(self):
        """Return m 1D arrays of length N, where m is the number of species
        in the Nanoparticle and N is the number of atoms.
        
        Each species has its own 1D array: V^E where E is the element of the specie.

        v^E_i = 1 if self[i].symbol == E 
        v^E_i = 0 if self[i].symbol != E
        
        Returns: np.arrays: 1D arrays that contains the occupancy of a lattice position
        base on the element of the array
        """
        self._occupancy_symbol_a = np.zeros(len(self))
        self._occupancy_symbol_b = np.zeros(len(self))

        for atom in self:
            if atom.symbol == list(self.symbols.indices())[0]:
                self._occupancy_symbol_a[atom.index] = 1
            else:
                self._occupancy_symbol_b[atom.index] = 1

        return self._occupancy_symbol_a, self._occupancy_symbol_b
    
    def get_connectivity_matrix(self):
        """Calculates the connectivity matrix A, an NxN matrix where N is the number of atoms.
        
            A_ij = 0 if atom_i and atom_j are not bonded together
            A_ij = 1 if atom_i and atom_j are bonded together

            Returns:
            np.array: Symmetric 2D matrix containing the connectivity between atoms.
        """
        from ase.neighborlist import natural_cutoffs
        from ase.neighborlist import build_neighbor_list

        self._connectivity_matrix = np.zeros((len(self), len(self)), dtype=np.int16)

        cutoffs = natural_cutoffs(self)    
        neighbor_list = build_neighbor_list(self,
                                            cutoffs=cutoffs,
                                            bothways=True,
                                            self_interaction=False)

        for atom_idx, _ in enumerate(self):
            neighbors, _ = neighbor_list.get_neighbors(atom_idx)
            for neighbor_idx in neighbors:
                self._connectivity_matrix[atom_idx][neighbor_idx] = 1

        return self._connectivity_matrix



