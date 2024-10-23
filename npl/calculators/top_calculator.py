
from typing import Union

from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.calculator import Calculator
import pickle

class TOPCalculator(Calculator):
    """
    A class representing a calculator for performing relaxation calculations using the ASE library.
    
    Parameters:
        calculator (Calculator): The calculator object used for performing the calculations.
        fmax (float): The maximum force tolerance for the relaxation.
    """
    
    def __init__(self,
                 model_paths : Union[list, str],
                 **kwargs
                 ):
        Calculator.__init__(self, **kwargs)
    
        self.model = pass
    
    def calculate(self, atoms):
        return NotImplementedError
