"""Definition of the canonical ensemble class."""

from ase import Atoms
from ase.units import kB
from typing import List

from .. import DataContainer
from ..calculators.base_calculator import BaseCalculator
from .thermodynamic_base_ensemble import ThermodynamicBaseEnsemble


class CanonicalEnsemble(ThermodynamicBaseEnsemble):
    """Instances of this class allow one to simulate systems in the
    canonical ensemble (:math:`N_iVT`), i.e. at constant temperature
    (:math:`T`), number of atoms of each species (:math:`N_i`), and
    volume (:math:`V`).

    The probability for a particular state in the canonical ensemble is
    proportional to the well-known Boltzmann factor,

    .. math::

        \\rho_{\\text{C}} \\propto \\exp [ - E / k_B T ].

    Since the concentrations or equivalently the number of atoms of each
    species is held fixed in the canonical ensemble, a trial step must
    conserve the concentrations. This is accomplished by randomly picking two
    unlike atoms and swapping their identities. The swap is accepted with
    probability

    .. math::

        P = \\min \\{ 1, \\, \\exp [ - \\Delta E / k_B T  ] \\},

    where :math:`\\Delta E` is the change in potential energy caused by the
    swap.

    The canonical ensemble provides an ideal framework for studying the
    properties of a system at a specific concentration. Properties such as
    potential energy or phenomena such as chemical ordering at a specific
    temperature can conveniently be studied by simulating at that temperature.
    The canonical ensemble is also a convenient tool for "optimizing" a
    system, i.e., finding its lowest energy chemical ordering. In practice,
    this is usually achieved by simulated annealing, i.e. the system is
    equilibrated at a high temperature, after which the temperature is
    continuously lowered until the acceptance probability is almost zero. In a
    well-behaved system, the chemical ordering at that point corresponds to a
    low-energy structure, possibly the global minimum at that particular
    concentration.

    Parameters
    ----------
    structure : :class:`Atoms <ase.Atoms>`
        atomic configuration to be used in the Monte Carlo simulation;
        also defines the initial occupation vector
    calculator : :class:`BaseCalculator <mchammer.calculators.ClusterExpansionCalculator>`
        calculator to be used for calculating the potential changes
        that enter the evaluation of the Metropolis criterion
    temperature : float
        temperature :math:`T` in appropriate units [commonly Kelvin]
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
        such as temperature or the number of atoms of different species
    trajectory_write_interval : int
        interval at which the current occupation vector of the atomic
        configuration is written to the data container.
    sublattice_probabilities : List[float]
        probability for picking a sublattice when doing a random swap.
        This should be as long as the number of sublattices and should
        sum up to 1.
    """

    def __init__(self,
                 structure: Atoms,
                 calculator: BaseCalculator,
                 temperature: float,
                 boltzmann_constant: float = kB,
                 random_seed: int = None,) -> None:

        self._ensemble_parameters = dict(temperature=temperature)

        # add species count to ensemble parameters
        symbols = set([symbol for sub in calculator.sublattices
                       for symbol in sub.chemical_symbols])
        for symbol in symbols:
            key = 'n_atoms_{}'.format(symbol)
            count = structure.get_chemical_symbols().count(symbol)
            self._ensemble_parameters[key] = count

        super().__init__(
            structure=structure,
            calculator=calculator,
            random_seed=random_seed,
            boltzmann_constant=boltzmann_constant)

    @property
    def temperature(self) -> float:
        """ Current temperature """
        return self._ensemble_parameters['temperature']

    def _do_trial_step(self):
        """ Carries out one Monte Carlo trial step. """
        
        sublattice_index = self.get_random_sublattice_index(self._swap_sublattice_probabilities)
        return self.do_canonical_swap(sublattice_index=sublattice_index)

