from abc import abstractproperty
from typing import Dict, List, Type, Any
import numpy as np

from ase import Atoms
from ase.units import kB
from ase.data import chemical_symbols

from npl.calculators.base_calculator import BaseCalculator
from npl.global_optimization.operations import BaseOperator
from npl.global_optimization.ensambles.base_ensamble import BaseEnsemble


class ThermodynamicBaseEnsemble(BaseEnsemble):
    """
    Thermodynamic base ensemble class.

    Parameters
    ----------
    structure : :class:`Atoms <ase.Atoms>`
        atomic configuration to be used in the Monte Carlo simulation;
        also defines the initial occupation vector
    calculator : :class:`BaseCalculator <mchammer.calculators.ClusterExpansionCalculator>`
        calculator to be used for calculating the potential changes
        that enter the evaluation of the Metropolis criterion
    boltzmann_constant : float
        Boltzmann constant :math:`k_B` in appropriate
        units, i.e. units that are consistent
        with the underlying cluster expansion
        and the temperature units [default: eV/K]
    user_tag : str
        human-readable tag for ensemble [default: None]
    random_seed : int
        seed for the random number generator used in the Monte Carlo
        simulation
    dc_filename : str
        name of file the data container associated with the ensemble
        will be written to; if the file exists it will be read, the
        data container will be appended, and the file will be
        updated/overwritten
    data_container_write_period : float
        period in units of seconds at which the data container is
        written to file; writing periodically to file provides both
        a way to examine the progress of the simulation and to back up
        the data [default: 600 s]
    ensemble_data_write_interval : int
        interval at which data is written to the data container; this
        includes for example the current value of the calculator
        (i.e. usually the energy) as well as ensembles specific fields
        such as temperature or the number of structure of different species
    trajectory_write_interval : int
        interval at which the current occupation vector of the atomic
        configuration is written to the data container.
    """

    def __init__(self,
                 structure: Atoms,
                 calculator: BaseCalculator,
                 operator: BaseOperator,
                 boltzmann_constant: float = kB,
                 random_seed: int = None,) -> None:

        self._boltzmann_constant = boltzmann_constant

        super().__init__(
            structure=structure,
            calculator=calculator,
            operator=operator,
            random_seed=random_seed,)

    @abstractproperty
    @property
    def temperature(self) -> float:
        pass

    @property
    def boltzmann_constant(self) -> float:
        """ Boltzmann constant :math:`k_B` (see parameters section above) """
        return self._boltzmann_constant

    def _acceptance_condition(self, potential_diff: float) -> bool:
        """
        Evaluates Metropolis acceptance criterion.

        Parameters
        ----------
        potential_diff
            change in the thermodynamic potential associated
            with the trial step
        """
        if potential_diff <= 0:
            return True
        elif self.temperature <= 1e-16:
            return False
        else:
            p = np.exp(-potential_diff / (self.boltzmann_constant * self.temperature))
            return p > self._next_random_number()

    def do_canonical_swap(self) -> int:
        """ Carries out one Monte Carlo trial step.

        Parameters
        ---------
        sublattice_index
            the sublattice the swap will act on
        allowed_species
            list of atomic numbers for allowed species

        Returns
        -------
        Returns 1 or 0 depending on if trial move was accepted or rejected
        """
        self.operator.perform_operation(self.structure)
        potential_diff = self._get_property_change()

        if self._acceptance_condition(potential_diff):
            return 1
        self.operator.revert_operation(self.structure)
        return 0

