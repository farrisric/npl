from typing import Union
import logging
from ase.calculators.calculator import Calculator
from sklearn.linear_model import LinearRegression
import npl.descriptors
import pickle
import numpy as np
import npl
import npl.descriptors

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
                 feature_classifier : npl.descriptors = None,
                 **kwargs
                 ):
        Calculator.__init__(self, **kwargs)

        self.feature_key = feature_key
        self.energy_key = feature_key

        if feature_classifier:
            self.feature_classifier = feature_classifier

        if model_paths:
            self.model = self.load_model(model_paths)

        if stoichiometry:
            self.coefficients = self.load_coefficients(stoichiometry)
            self.model = LinearRegression()
            self.model.coef_ = self.coefficients

    def load_coefficients(self, stoichiometry):
        logging.info("Loading top parameters of {}".format(stoichiometry))

        params = self.get_data_by_stoichiometry(stoichiometry)
        self.feature_classifier = self.feature_classifier(params['symbols'])
        feature_name = self.feature_classifier.get_feature_labels()
        coefficients = np.zeros(len(feature_name))

        for i, feature in enumerate(feature_name):
            coefficients[i] = params.get(feature, 0)

        for i, feature in enumerate(feature_name):
            if 'cn' in feature:
                symbol = feature[:2]
                cn = feature[3:-1]
                coefficients[i] = params['data'][symbol][cn]

        logging.info("Parameters obtained from reference: {}".format(params['reference']))
        logging.info("Parameters loaded successfully")
        non_zero_params = {feature_name[i]: coef for i, coef in enumerate(coefficients)
                           if coef != 0}
        logging.info("Parameters: \n{}".format(non_zero_params))
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
        from .parameters import top_parameters
        return top_parameters[stoichiometry]

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

    def get_feature_classifier(self):
        return self.feature_classifier
