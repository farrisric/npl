from ..mc_moves import InsertionMove, DeletionMove, DisplacementMove
from .canonical_ensemble import BaseEnsemble
from npl.utils.random_number_generator import RandomNumberGenerator
import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator
from typing import Optional
import logging


PLANCK_CONSTANT = 1  # 6.62607e-34  # Planck's constant in m²kg/s
BOLTZMANN_CONSTANT = 1  # 1.38066e-23  # Boltzmann constant in J/K

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class GrandCanonicalEnsemble(BaseEnsemble):
    def __init__(self,
                 atoms: Atoms,
                 calculator: Calculator,
                 mu : dict,
                 masses: dict,
                 species : list,
                 temperature : float,
                 moves: dict,
                 user_tag: Optional[str] = None,
                 random_seed: Optional[int] = None,
                 trajectory_write_interval: Optional[int] = None,
                 traj_file: str = 'traj_test.traj',
                 outfile: str = 'outfile.out',
                 outfile_write_interval: int = 10) -> None:
        
        super().__init__(atoms,
                         calculator,
                         random_seed,
                         traj_file,
                         outfile,
                         outfile_write_interval)

        self.volume = atoms.get_volume()
        self.masses = masses
        self.n_atoms = len(self.atoms)
        self.species = species
        self._temperature = temperature
        self._mu = mu
        self._beta = 1/(self._temperature*BOLTZMANN_CONSTANT)

        self.n_ins_del = moves[0]
        self.n_displ = moves[1]
        self.n_moves = self.n_ins_del + self.n_displ

        self.frac_ins_del = self.n_ins_del/self.n_moves
        self.rng_move_choice = RandomNumberGenerator(seed=self._random_seed)
        self.rng_acceptance = RandomNumberGenerator(seed=self._random_seed)
        # Initialize GCMC moves
        self.insert_move = InsertionMove(species=self.species, seed=self._random_seed)
        self.deletion_move = DeletionMove(species=self.species, seed=self._random_seed)
        self.displace_move = DisplacementMove(species=self.species, seed=self._random_seed)

        # COUNTERS
        self.count_moves = {'Displacements' : 0, 'Insertions' : 0, 'Deletions' : 0}
        self.count_acceptance = {'Displacements' : 0, 'Insertions' : 0, 'Deletions' : 0}

    def _acceptance_condition(self,
                              potential_diff: float,
                              delta_particles: int,
                              species: str) -> bool:
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
            if potential_diff <= 0:
                return True
            else:
                p = np.exp(-potential_diff * self._beta)
                return p > self.rng_acceptance.get_uniform()

        # de Broglie wavelength in meters
        lambda_db = 1  # PLANCK_CONSTANT / np.sqrt(2 * np.pi * self.masses[species] * (1 / self._beta))

        if delta_particles == 1:  # Insertion move
            db_term = (self.volume / ((self.n_atoms+1)*lambda_db**3))
            p = db_term * np.exp(-self._beta * (potential_diff - self._mu[species]))

        elif delta_particles == -1:  # Deletion move
            db_term = (lambda_db**3*self.n_atoms / self.volume)
            p = db_term * np.exp(-self._beta * (potential_diff + self._mu[species]))

        if p > 1:
            return True
        else:
            return p > self.rng_acceptance.get_uniform()

    def do_displ_move(self):
        self.count_moves['Displacements'] += 1
        return self.displace_move.do_trial_move(self.atoms)

    def do_ins_del_move(self):
        r2 = self.rng_move_choice.get_uniform()
        if r2 <= 0.5:
            self.count_moves['Insertions'] += 1
            return self.insert_move.do_trial_move(self.atoms)
        else:
            self.count_moves['Deletions'] += 1
            return self.deletion_move.do_trial_move(self.atoms)

    def do_trial_step(self, ):
        for move in range(self.n_moves):
            r1 = self.rng_move_choice.get_uniform()
            if r1 <= self.frac_ins_del:
                atoms_new, delta_particles, species = self.do_ins_del_move()
            else:
                atoms_new = self.do_displ_move()
                delta_particles = 0
                species = 'X'

            E_new = self.compute_energy(atoms_new)
            delta_E = E_new - self.E_old
            if self._acceptance_condition(delta_E, delta_particles, species):
                self.atoms = atoms_new
                self.n_atoms = len(self.atoms)
                self.E_old = E_new
                if delta_particles == 0:
                    self.count_acceptance['Displacements'] += 1
                if delta_particles == 1:
                    self.count_acceptance['Insertions'] += 1
                if delta_particles == -1:
                    self.count_acceptance['Deletions'] += 1

    def compute_energy(self, atoms):
        return self._calculator.get_potential_energy(atoms)

    def run(self, steps):
        """
        Runs the Grand Canonical Monte Carlo simulation.

        Returns:
            None
        """
        logger.info('+---------------------------------+')
        logger.info('| Grand Canonical Ensemble Monte Carlo  |')
        logger.info('+---------------------------------+')
        logger.info('Starting simulation...')
        logger.info('Temperature: {}'.format(self._temperature))
        logger.info('Number of steps: {}'.format(steps))

        self.E_old = self.compute_energy(self.atoms)
        #self.write_traj_file(self.atoms)
        #self.write_outfile(self._step, self.lowest_energy)

        for _ in range(steps):
            accepted = self.do_trial_step()
            
            self._step += 1
            # self._accepted_trials += accepted

            if self._step % self._outfile_write_interval == 0:
                print(self._step,
                      self.n_atoms,
                      self.E_old,
                      np.array(list(self.count_acceptance.values()))/np.array(list(self.count_moves.values())))
            #     self.write_outfile(self._step, self.lowest_energy)
            #     logger.info('Step: {}'.format(self._step))
            #     logger.info('Lowest energy: {}'.format(self.lowest_energy))
            #     logger.info('Accepted trials: {}'.format(self._accepted_trials))
