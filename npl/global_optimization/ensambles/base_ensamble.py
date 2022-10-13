from abc import ABC, abstractmethod

from npl.calculators import BaseCalculator

from ase import Atoms


class BaseEnsemble(ABC):
    """Base ensemble class.

    Parameters
    ----------
    structure : :class:`Atoms <ase.Atoms>`
        atomic configuration to be used in the Monte Carlo simulation;
        also defines the initial occupation vector
    calculator : :class:`BaseCalculator <mchammer.calculators.ClusterExpansionCalculator>`
        calculator to be used for calculating the potential changes
        that enter the evaluation of the Metropolis criterion
    random_seed : int
        seed for the random number generator used in the Monte Carlo
        simulation
    """

    def __init__(self,
                 structure: Atoms,
                 calculator: BaseCalculator,
                 random_seed: int = None,) -> None:

        # initialize basic variables
        self._accepted_trials = 0
        self._observers = {}  # type: Dict[str, BaseObserver]
        self._step = 0

        # calculator and configuration
        self._calculator = calculator
