import random

from abc import ABC, abstractmethod

from ase import Atoms

from npl.calculators.base_calculator import BaseCalculator
from npl.global_optimization.operations import BaseOperator

class BaseEnsemble(ABC):
    """Base ensemble class.

    Parameters
    ----------
    structure : :class:`Atoms <ase.Atoms>`
        atomic configuration to be used in the Monte Carlo simulation;
        also defines the initial occupation vector
    calculator : :class:`BaseCalculator <npl.calculators.TopologicalCalculator>`
        calculator to be used for calculating the potential changes
        that enter the evaluation of the Metropolis criterion
    random_seed : int
        seed for the random number generator used in the Monte Carlo
        simulation
    operaton : :class: `BaseOperation <npl.operations.BaseOperator>`
        operation applied to every Monte Carlo step;
    """

    def __init__(self,
                 structure: Atoms,
                 calculator: BaseCalculator,
                 operator: BaseOperator,
                 random_seed: int = None) -> None:

        # initialize basic variables
        self._accepted_trials = 0
        self._observers = {}  # type: Dict[str, BaseObserver]
        self._step = 0

        # calculator and configuration
        self._calculator = calculator

        if random_seed is None:
            self._random_seed = random.randint(0, int(1e16))
        else:
            self._random_seed = random_seed
        random.seed(a=self._random_seed)


    @property
    def structure(self) -> Atoms:
        """ current configuration (copy) """
        return self.configuration.structure

    @property
    def calculator(self) -> BaseCalculator:
        """ calculator attached to the ensemble """
        return self._calculator

    @property
    def step(self) -> int:
        """ current trial step counter """
        return self._step


    def _run(self, number_of_trial_steps: int):
        """Runs MC simulation for a number of trial steps without
        interruption.

        Parameters
        ----------
        number_of_trial_steps
            number of trial steps to run without stopping
        """
        for _ in range(number_of_trial_steps):
            accepted = self._do_trial_step()
            self._step += 1
            self._accepted_trials += accepted

    @abstractmethod
    def _do_trial_step(self):
        pass




