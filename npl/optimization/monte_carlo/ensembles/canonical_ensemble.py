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
                 traj_file='traj_test.traj', 
                 op_list=None, 
                 outfile='outfile.out',
                 outfile_write_interval=10,
                 constraints=None) -> None:

        super().__init__(structure=atoms, calculator=calculator, random_seed=random_seed)

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
        self._traj_file = traj_file
        self._op_list = op_list
        self._outfile = outfile
        self._outfile_write_interval = outfile_write_interval

        self._step = 0
        self._accepted_trials = 0

    def write_outfile(self, step, energy):
        with open(self._outfile, 'a') as outfile:
            outfile.write(' STEP: {} ENERGY: {}\n'.format(step, energy))
    
    def write_traj_file(self, atoms):
        Trajectory(self._traj_file, 'a').write(atoms)

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


if __name__ == '__main__':
    from ase.io import read
    from ase.optimize import BFGS
    from mace.calculators import mace_mp
    from acat.ga.adsorbate_operators import (AddAdsorbate, RemoveAdsorbate,
                                         MoveAdsorbate, ReplaceAdsorbate,
                                         SimpleCutSpliceCrossoverWithAdsorbates)
    from acat.ga.particle_mutations import (RandomPermutation, COM2surfPermutation,
                                            Rich2poorPermutation, Poor2richPermutation)
    from ase.ga.offspring_creator import OperationSelector
    from acat.adsorption_sites import ClusterAdsorptionSites
    from acat.adsorbate_coverage import ClusterAdsorbateCoverage

    atoms = read('/home/g15farris/AdsGO/PdAg_CO_H/bare_nps/xyz/Ag300Pd105_LEH.xyz')
    atoms.center(5)
    calculator = mace_mp(model="/work/g15farris/2023-12-03-mace-128-L1_epoch-199.model")
    optimizer = BFGS
    temperature = 300
    steps = 100
    traj_file = 'traj_test.traj'
    sas = ClusterAdsorptionSites(atoms, composition_effect=True)
    soclist = ([1, 1, 1, 1, 1],
        [Rich2poorPermutation(elements=['Ag', 'Pd'], num_muts=1),
         Poor2richPermutation(elements=['Ag', 'Pd'], num_muts=1),
         RandomPermutation(elements=['Ag', 'Pd'], num_muts=1),
         COM2surfPermutation(elements=['Ag', 'Pd'], num_muts=1),
         MoveAdsorbate(['H'], adsorption_sites=sas, num_muts=1)
        ])
    op_list = OperationSelector(*soclist)

    outfile = 'outfile.out'
    outfile_write_interval = 10

    montecarlo = CanonicalEnsemble(atoms,
                 calculator,
                 optimizer, 
                 fmax=0.1,
                 temperature=temperature, 
                 steps=steps, 
                 traj_file=traj_file, 
                 op_list=op_list, 
                 outfile=outfile,
                 outfile_write_interval=outfile_write_interval)
    
    montecarlo.run()
