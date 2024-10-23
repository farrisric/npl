from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.calculator import Calculator
from base_calculator import BaseCalculator

class Calculator(BaseCalculator):
    """
    A class representing a calculator for performing relaxation calculations using the ASE library.
    
    Parameters:
        calculator (Calculator): The calculator object used for performing the calculations.
        fmax (float): The maximum force tolerance for the relaxation.
    """
    
    def __init__(self,
                 calculator : Calculator,
                 fmax : float,
                 ):
        super().__init__()
        self._calculator = calculator
        self._fmax = fmax

    def relax(self, structure : Atoms):
        """
        Perform relaxation of the given structure using the specified calculator and force tolerance.
        
        Parameters:
            structure (Atoms): The structure to be relaxed.
        
        Returns:
            tuple: A tuple containing the relaxed structure and the potential energy.
        """
        structure.calc = self._calculator
        optimizer = BFGS(structure)
        optimizer.run(fmax=self._fmax)
        
        return structure, structure.get_potential_energy()
