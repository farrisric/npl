import random

from abc import ABC, abstractmethod
from time import time

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

    def _terminate_sampling(self) -> bool:
        """This method is called from the run method to determine whether the MC
        sampling loop should be terminated for a reason other than having exhausted
        the number of iterations. The method can be overriden by child classes in
        order to provide an alternative exit mechanism.
        """
        return False

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
        
    @property
    def random_seed(self) -> int:
        """ seed used to initialize random number generator """
        return self._random_seed

    def _next_random_number(self) -> float:
        """ Returns the next random number from the PRNG. """
        return random.random()


    def run(self, number_of_trial_steps: int):
        """
        Samples the ensemble for the given number of trial steps.

        Parameters
        ----------
        number_of_trial_steps
            number of MC trial steps to run in total

        Raises
        ------
        TypeError
            if `number_of_trial_steps` is not an int
        """

        if not isinstance(number_of_trial_steps, int):
            raise TypeError('number_of_trial_steps must be an integer ({})'
                            .format(number_of_trial_steps))

        last_write_time = time()

        initial_step = self.step
        final_step = self.step + number_of_trial_steps
        # run Monte Carlo simulation such that we start at an
        # interval which lands on the observer interval
        if initial_step != 0:
            first_run_interval = self.observer_interval -\
                (initial_step -
                 (initial_step // self.observer_interval) *
                 self.observer_interval)
            first_run_interval = min(first_run_interval, number_of_trial_steps)
            self._run(first_run_interval)
            initial_step += first_run_interval

        step = initial_step
        while step < final_step and not self._terminate_sampling():
            uninterrupted_steps = min(self.observer_interval, final_step - step)
            if self.step % self.observer_interval == 0:
                self._observe(self.step)

            self._run(uninterrupted_steps)
            step += uninterrupted_steps

        # if we end on an observation interval we also observe
        if self.step % self.observer_interval == 0:
            self._observe(self.step)

        # allow ensemble a chance to go clean
        #self._finalize()


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

    def _get_property_change(self) -> float:
        """Computes and returns the property change due to a change of
        the configuration.

        _N.B.:_ This method leaves the configuration itself unchanged.

        Parameters
        ----------
        sites
            indices of sites to change
        species
            new occupations (species) by atomic number
        """
        return self.calculator.calculate_change()





