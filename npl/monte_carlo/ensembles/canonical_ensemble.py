from ase.units import kB as boltzmann_constant
from ase import Atoms
from ase.io.trajectory import Trajectory
import numpy as np
import random
from .base_ensemble import BaseEnsemble

class CanonicalEnsemble(BaseEnsemble):
    def __init__(self,
                 atoms,
                 calculator,
                 random_seed=None,
                 optimizer=None,
                 fmax=0.1,
                 temperature=300, 
                 steps=100, 
                 op_list=None, 
                 constraints=None,
                 traj_file: str = 'traj_test.traj', 
                 outfile: str = 'outfile.out',
                 outfile_write_interval: int = 10) -> None:

        super().__init__(structure=atoms, calculator=calculator, random_seed=random_seed, traj_file=traj_file,outfile=outfile,outfile_write_interval=outfile_write_interval)

        if random_seed is not None:
            random.seed(random_seed)

        self.lowest_energy = float('inf')  # Initialize to positive infinity
        self.atoms = atoms
        self.constraints = constraints
        self._temperature = temperature
        self._calculator = calculator
        self._optimizer = optimizer
        self._fmax = fmax
        self._steps = steps
        self._op_list = op_list

        self._step = 0
        self._accepted_trials = 0

    def _acceptance_condition(self, potential_diff: float) -> bool:
        if potential_diff <= 0:
            return True
        elif self._temperature <= 1e-16:
            return False
        else:
            p = np.exp(-potential_diff / (boltzmann_constant * self._temperature))
            return p > self._next_random_number()
    
    def _next_random_number(self) -> float:
        """Returns the next random number from the PRNG."""
        return random.random()

    def relax(self, atoms) -> Atoms:
        atoms.info['key_value_pairs'] = {}
        atoms.calc = self._calculator
        opt = self._optimizer(atoms, logfile=None)
        opt.run(fmax=self._fmax)
        
        Epot = atoms.get_potential_energy()
        atoms.info['key_value_pairs']['potential_energy'] = Epot
        return atoms

    def do_mutations(self):
        new_atoms = self.atoms.copy()  # Use copy to avoid modifying the original
        for op_i in range(np.random.geometric(p=0.35, size=1)[0]):
            new_atoms.info['data'] = {'tag': None}
            new_atoms.info['confid'] = 1
            operation = self._op_list.get_operator()
            new_atoms, _ = operation.get_new_individual([new_atoms])
        if self.constraints:
            new_atoms.set_constraint(self.constraints)    
        new_atoms = self.relax(new_atoms)
        return new_atoms

    def trial_step(self):
        new_atoms = self.do_mutations()
        potential_diff = new_atoms.info['key_value_pairs']['potential_energy'] - self.atoms.info['key_value_pairs']['potential_energy']
        
        if self._acceptance_condition(potential_diff):
            if new_atoms.info['key_value_pairs']['potential_energy'] < self.lowest_energy:
                self.lowest_energy = new_atoms.info['key_value_pairs']['potential_energy']
            self.atoms = new_atoms
            self.write_traj_file(self.atoms)
            return 1
        return 0

    def run(self):
        self.relax(self.atoms)
        self.lowest_energy = self.atoms.get_potential_energy()
        for _ in range(self._steps):
            accepted = self.trial_step()
            self._step += 1
            self._accepted_trials += accepted

            if self._step % self._outfile_write_interval == 0:
                self.write_outfile(self._step, self.lowest_energy)

