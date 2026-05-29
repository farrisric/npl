import copy
import pickle


class EnergyCalculator:
    """Base class for an energy calculator.

    Valid implementations have to implement the compute_energy(particle) function.
    Energies are saved in the particle object with the key of the respective calculator.
    """
    def __init__(self):
        self.energy_key = None
        pass

    def compute_energy(self, particle):
        raise NotImplementedError

    def get_energy_key(self):
        return copy.deepcopy(self.energy_key)

    def set_energy_key(self, energy_key):
        self.energy_key = energy_key

    def save(self, name_file : str):
        with open(name_file, 'wb') as out:
            pickle.dump(self, out)

    @staticmethod
    def load(name_file):
        with open(name_file, 'rb') as calc:
            return pickle.load(calc)
