from abc import ABC, abstractmethod
from ase import Atoms
from npl.utils import RandomNumberGenerator
import numpy as np


class BaseMove(ABC):
    """Abstract base class for Monte Carlo moves."""

    def __init__(self, species: list[str], rng: RandomNumberGenerator) -> None:
        """
        Initializes the move with the given atomic configuration, species, and RNG.

        Parameters:
        atoms (Atoms): ASE Atoms object representing the system.
        species (list[str]): List of possible atomic species for insertion.
        rng (RandomNumberGenerator): Random number generator.
        """
        self.species = species
        self.rng = rng

    @abstractmethod
    def do_trial_move(self, atoms) -> Atoms:
        """
        Perform the Monte Carlo move and return the new atomic configuration.

        Returns:
        Atoms: Updated ASE Atoms object after the move.
        """
        pass


class InsertionMove(BaseMove):
    """Class for performing an insertion move."""

    def do_trial_move(self, atoms) -> Atoms:
        """
        Insert a random atom of a random species at a random position.

        Returns:
        Atoms: Updated ASE Atoms object after the insertion.
        """
        box = atoms.get_cell()
        selected_species = self.rng.random.choice(self.species)
        position = np.array([
            box[i]*self.rng.get_uniform() for i in range(3)
            ]).sum(axis=1)
        return Atoms(selected_species, positions=[position])


class DeletionMove(BaseMove):
    """Class for performing a deletion move."""

    def do_trial_move(self, atoms) -> int:
        """
        Delete a random atom from the structure.

        Returns:
        Atoms: Updated ASE Atoms object after the deletion.
        """
        selected_species = self.rng.random.choice(self.species)
        indices_of_species = [atom.index for atom in atoms if atom.symbol in selected_species]
        if len(indices_of_species) == 0:
            return False
        atom_index = self.rng.random.choice(indices_of_species)
        return atom_index


class DisplacementMove(BaseMove):
    """Class for performing a displacement move."""

    def __init__(self,
                 species: list[str],
                 rng: RandomNumberGenerator, max_displacement: float = 0.1) -> None:
        """
        Initializes the displacement move with a maximum displacement.

        Parameters:
        max_displacement (float): Maximum displacement distance.
        """
        super().__init__(species, rng)
        self.max_displacement = max_displacement

    def do_trial_move(self, atoms) -> Atoms:
        """
        Displace a random atom by a random vector within the maximum displacement range.

        Returns:
        Atoms: Updated ASE Atoms object after the displacement.
        """
        if len(atoms) == 0:
            raise ValueError("No atoms to displace.")
        atom_index = self.rng.random.randint(0, len(atoms) - 1)
        displacement = [
            self.rng.get_uniform(-self.max_displacement, self.max_displacement) for _ in range(3)
            ]
        return atom_index, displacement
