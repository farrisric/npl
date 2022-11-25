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

        self._connectivity_matrix = self.get_connectivity_matrix()
        self._occupation_matrix = self.get_occupation_matrix()
        self.neighbor_dict = {x.index : [] for x in self}

        self.constructu_neighbor_list()

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
    
    def get_bond_matrix(self):
        return np.dot(self._connectivity_matrix, self._occupation_matrix.T)

    def get_occupation_matrix(self):
        """Return m NxM matrix where M is the number of species
        in the Nanoparticle and N is the number of atoms.
        
        Each species has its own 1D array: V^E where E is the element of the specie.

        v^E_i = 1 if self[i].symbol == E 
        v^E_i = 0 if self[i].symbol != E
        
        Returns: np.arrays: 1D arrays that contains the occupancy of a lattice position
        base on the element of the array
        """
        elements = sorted(list(self.symbols.indices()))
        n_elements = len(elements)

        occupation_matrix = np.zeros((n_elements, self.get_global_number_of_atoms()))

        for atom in self:
            element_index = elements.index(atom.symbol)
            occupation_matrix[element_index, atom.index] = 1

        return occupation_matrix
    
    def get_connectivity_matrix(self):
        """Calculates the connectivity matrix A, an NxN matrix where N is the number of atoms.
        
            A_ij = 0 if atom_i and atom_j are not bonded together
            A_ij = 1 if atom_i and atom_j are bonded together

            Returns:
            np.array: Symmetric 2D matrix containing the connectivity between atoms.
        """
        from ase.neighborlist import natural_cutoffs
        from ase.neighborlist import build_neighbor_list

        connectivity_matrix = np.zeros((len(self), len(self)), dtype=np.int16)

        cutoffs = natural_cutoffs(self)    
        neighbor_list = build_neighbor_list(self,
                                            cutoffs=cutoffs,
                                            bothways=True,
                                            self_interaction=False)

        for atom_idx, _ in enumerate(self):
            neighbors, _ = neighbor_list.get_neighbors(atom_idx)
            for neighbor_idx in neighbors:
                connectivity_matrix[atom_idx][neighbor_idx] = 1

        return connectivity_matrix

    def constructu_neighbor_list(self):
        """Return a dictionary where the keys are the atom indices and the
        values are the indices of its first nearest-neighbors."""
        from ase.neighborlist import natural_cutoffs
        from ase.neighborlist import build_neighbor_list

        self.neighbor_dict = {x.index : [] for x in self}
        cutoffs = natural_cutoffs(self)    
        neighbor_list = build_neighbor_list(self,
                                            cutoffs=cutoffs,
                                            bothways=True,
                                            self_interaction=False)

        for atom_idx in self.neighbor_dict:
            self.neighbor_dict[atom_idx] = neighbor_list.get_neighbors(atom_idx)[0]

    def get_coordination_list(self):
        coordination_list = np.empty(len(self))

        for index, neigh in self.neighbor_dict.items():
            coordination_list[index] = len(neigh)

        return coordination_list

    def get_surface_neighbors(self):
        surface_indices = np.where(self.get_coordination_list()<12)[0]

        surface_neighbors_dict = {x : [] for x in surface_indices}
        
        for surface_idx in surface_indices:
            for neigh in self.neighbor_dict[surface_idx]:
                if neigh in surface_indices:
                    surface_neighbors_dict[surface_idx].append(neigh)

        return surface_neighbors_dict

        




        




