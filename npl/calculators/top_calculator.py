from typing import Union
import logging
from ase.calculators.calculator import Calculator
from sklearn.linear_model import LinearRegression
from npl.descriptors import TopologicalFeatureClassifier as TOP
import pickle
import json
import numpy as np
import os
import npl

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class TOPCalculator(Calculator):
    """
    A class representing a calculator for performing relaxation calculations using the ASE library.

    Parameters:
        calculator (Calculator): The calculator object used for performing the calculations.
        fmax (float): The maximum force tolerance for the relaxation.
    """

    def __init__(self,
                 feature_key : str,
                 stoichiometry : str = None,
                 model_paths : Union[list, str] = None,
                 **kwargs
                 ):
        Calculator.__init__(self, **kwargs)

        self.feature_key = feature_key
        self.energy_key = feature_key

        if model_paths:
            self.model = self.load_model(model_paths)

        if stoichiometry:
            self.coefficients = self.load_coefficients(stoichiometry)
            self.model = LinearRegression()
            self.model.coef_ = self.coefficients

    def load_coefficients(self, stoichiometry):
        logging.info("Loading top parameters of {}".format(stoichiometry))

        params = self.get_data_by_stoichiometry(stoichiometry)
        top = TOP(params['symbols'])
        feature_name = top.get_feature_labels()
        coefficients = np.zeros(len(feature_name))

        for i, feature in enumerate(feature_name):
            coefficients[i] = params.get(feature, 0)

        for i, feature in enumerate(feature_name):
            if len(feature) > 4:
                symbol = feature[:2]
                cn = feature[3:-1]
                coefficients[i] = params['data'][symbol][cn]

        logging.info("Parameters obtained from reference: {}".format(params['reference']))
        return coefficients

    def get_data_by_stoichiometry(self, stoichiometry):
        """
        Retrieve data for a given stoichiometry.

        Parameters:
            stoichiometry (str): The stoichiometry to look up (e.g., "Pt70Au70").
            data (dict): The JSON data.

        Returns:
            dict: The data for the specified stoichiometry, or None if not found.
        """
        params = os.path.join(os.path.dirname(npl.__file__), '../data', 'params.json')
        data = json.load(open(params))
        logging.debug("Data loaded from {}".format(params))
        return data.get(stoichiometry, None)

    def load_model(self, model_path):
        with open(model_path, 'rb') as calc:
            return pickle.load(calc)

    def calculate(self, atoms):
        self.results = {}
        feature_vector = atoms.info[self.feature_key]
        self.results['energy'] = self.model.predict(feature_vector)

    def compute_energy(self, particle):
        feature_vector = particle.get_feature_vector(self.feature_key)
        top_energy = np.dot(np.transpose(self.model.coef_), feature_vector)
        particle.set_energy(self.energy_key, top_energy)
        return top_energy

    def set_coefficients(self, coefficients):
        self.coefficients = coefficients
        self.model.coef_ = self.coefficients

    def get_coefficients(self):
        return self.coefficients

    def get_energy_key(self):
        return self.energy_key

    def get_feature_key(self):
        return self.feature_key
