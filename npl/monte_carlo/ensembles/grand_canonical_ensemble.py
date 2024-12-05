from .canonical_ensemble import CanonicalEnsemble
import numpy as np

PLANCK_CONSTANT = 6.62607e-34  # Planck's constant in m²kg/s
BOLTZMANN_CONSTANT = 1.38066e-23  # Boltzmann constant in J/K

class GrandCanonicalEnsemble(CanonicalEnsemble):
    def __init__(self,
                 atoms,
                 calculator,
                 random_seed=None,
                 optimizer=None,
                 fmax=0.1,
                 temperature=300,
                 op_list=None,
                 constraints=None,
                 mu,
                 p=1,
                 traj_file: str = 'traj_test.traj',
                 outfile: str = 'outfile.out',
                 outfile_write_interval: int = 10) -> None:
        super().__init__(atoms,
                         calculator,
                         random_seed,
                         optimizer,
                         fmax,
                         temperature,
                         op_list,
                         constraints,
                         p,
                         traj_file,
                         outfile,
                         outfile_write_interval)

        self.volume = atoms.get_volume()
        self.n_atoms = len(self.atoms)
        self.beta = 1/(self._temperature*BOLTZMANN_CONSTANT)
        self.mu = mu
        
    def _acceptance_condition(self, mass : float,
                              potential_diff: float, delta_particles: int,
                              mu : float) -> bool:
        """
        Determines whether to accept a trial move based on the potential energy difference and
        temperature.

        Args:
            potential_diff (float): The potential energy difference between the current and new
            configurations.

        Returns:
            bool: True if the trial move is accepted, False otherwise.
        """
        if delta_particles == 0:
            p = np.exp(-potential_diff / (PLANCK_CONSTANT * self._temperature))
            return p > self._next_random_number()

        lambda_db = PLANCK_CONSTANT / np.sqrt(2 * np.pi * mass * (1 / self.beta))  # de Broglie wavelength in meters

        if delta_particles == 1:  # Insertion move
            db_term = (self.volume / ((self.n_atoms+1)*lambda_db**3))
            probability = db_term * np.exp(-self.beta * (potential_diff - mu))
        elif delta_particles == -1:  # Deletion move
            db_term = (lambda_db**3*self.n_atoms / self.volume)
            probability = db_term * np.exp(-self.beta * (potential_diff + mu))
        return min(1.0, probability)
