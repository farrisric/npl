from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.calculator import Calculator
from base_calculator import BaseCalculator

class Calculator(BaseCalculator):
    
    def __init__(self,
                 calculator : Calculator,
                 fmax : float,
                 ):
        super().__init__()
        self._calculator = calculator
        self._fmax = fmax

    def relax(self, structure : Atoms):
        structure.calc = self._calculator
        optimizer = BFGS(structure)
        optimizer.run(fmax=self._fmax)
        
        return structure, structure.get_potential_energy()
