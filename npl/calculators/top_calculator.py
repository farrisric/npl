
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
                 feature_key : str,
                 **kwargs
                 ):
        Calculator.__init__(self, **kwargs)

        self.feature_key = feature_key
        self.model = self.load_model(model_paths)
    
    def load_model(self, model_path):
        with open(model_path, 'rb') as calc:
            return pickle.load(calc)
    
    def calculate(self, atoms):
        self.results = {}
        feature_vector = atoms.info[self.feature_key]
        self.results['energy'] = self.model.predict(feature_vector)
