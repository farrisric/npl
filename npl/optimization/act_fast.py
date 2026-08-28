"""Rust fast path for the ACT guided-exchange ordering search.

The ACT inner loop is integer bookkeeping plus one dot product, but in Python
it costs ~30 million interpreter calls per basin-hopping restart: profiling a
2000-attempt restart on a 405-atom particle put 19 s of 19 s inside
``ACTExchangeOperator.update`` and the local-environment recompute below it.
``npl_act_rs`` reimplements exactly that loop; this module packs a particle into
the flat arrays it wants and unpacks the result.

The fast path only takes over when every assumption behind it holds - two
elements, a linear TEC model, the extension actually importable - and
``run_basin_hopping``/``local_optimization`` fall back to the Python otherwise,
so nothing here changes what npl can do, only how fast it does it.
"""

import copy
import os

import numpy as np

from ase.data import atomic_numbers

try:  # the extension is optional; npl runs without it
    import npl_act_rs as _rs
except ImportError:  # pragma: no cover - exercised by installs without the wheel
    _rs = None


TEC_N_FEATURES = 182

#: Set this in the environment to keep a whole run on the pure Python path,
#: which is what a long unattended sweep wants if the extension lands
#: mid-flight: every process it spawns then uses the same code.
DISABLE_ENV_VAR = 'NPL_DISABLE_ACT_RS'


def available():
    """True when the compiled extension can be used."""
    return _rs is not None and not os.environ.get(DISABLE_ENV_VAR)


def usable(particle, energy_calculator, environment_energies, model):
    """Whether this call is one the Rust loop reproduces exactly.

    Every condition here is an assumption baked into the port: the ACT operator,
    two elements, and an energy that is a dot product of ridge coefficients with
    the 182-bin TEC histogram.
    """
    if not available() or model != 'ACT':
        return False
    if getattr(energy_calculator, 'feature_key', None) != 'TEC':
        return False
    if not hasattr(energy_calculator, 'get_coefficients'):
        return False
    if len(set(particle.get_all_symbols())) != 2:
        return False
    if environment_energies is None:
        return False
    if len(environment_energies) != TEC_N_FEATURES:
        return False
    try:
        coefficients = energy_calculator.get_coefficients()
    except Exception:
        return False
    return coefficients is not None and len(coefficients) == TEC_N_FEATURES


def pack(particle):
    """Flatten a particle into (offsets, neighbours, symbol codes, elements).

    Symbol code 0 is the alphabetically first element, matching
    ``TopologicalEnvironmentClassifier``, whose feature is
    ``91 * symbols.index(symbol) + cn_offset[cn] + n_first_element_neighbours``.
    Neighbour sets are emitted in sorted order; the features and exchange
    energies are sums over the set, so the order is free.
    """
    n_atoms = len(particle.get_indices())
    offsets = np.zeros(n_atoms + 1, dtype=np.uint32)
    neighbours = []
    for index in range(n_atoms):
        neighbour_indices = sorted(particle.neighbor_list[index])
        neighbours.extend(neighbour_indices)
        offsets[index + 1] = len(neighbours)

    elements = sorted(particle.get_all_symbols())
    first_number = atomic_numbers[elements[0]]
    numbers = np.asarray(particle.atoms.atoms.numbers)
    codes = np.where(numbers == first_number, 0, 1).astype(np.uint8)

    return (offsets.tolist(), [int(x) for x in neighbours], codes.tolist(),
            elements)


def unpack(particle, codes, elements):
    """Write a symbol-code vector back onto a copy of the particle."""
    numbers = np.array([atomic_numbers[elements[code]] for code in codes])
    particle.atoms.atoms.set_atomic_numbers(numbers)
    return particle


def _coefficients(energy_calculator):
    return [float(c) for c in energy_calculator.get_coefficients()]


def _restore(particle, local_env_calculator, local_feature_classifier):
    """Recompute environments and features so the particle is self-consistent
    with whatever ordering it now carries."""
    local_env_calculator.compute_local_environments(particle)
    local_feature_classifier.compute_feature_vector(particle)


def basin_hopping(start_particle, energy_calculator, environment_energies,
                  n_hopping_attempts, n_hops, local_env_calculator,
                  local_feature_classifier, seed):
    """Rust equivalent of ``run_basin_hopping`` for the ACT case.

    Returns the same ``[best_particle, lowest_energies, flip_energy_list]`` the
    Python returns.
    """
    offsets, neighbours, codes, elements = pack(start_particle)
    best_codes, best_energy, lowest_energies, flips, error = _rs.act_basin_hopping(
        offsets, neighbours, codes, _coefficients(energy_calculator),
        [float(e) for e in environment_energies],
        int(n_hopping_attempts), int(n_hops), int(seed))

    best_particle = unpack(start_particle, best_codes, elements)
    _restore(best_particle, local_env_calculator, local_feature_classifier)
    best_particle.set_energy(energy_calculator.get_energy_key(), best_energy)

    flip_energy_list = [(flip, delta, list(pair)) for flip, delta, pair in flips]
    return [best_particle, lowest_energies, flip_energy_list], error


def genetic_algorithm(start_population, energy_calculator, environment_energies,
                      local_env_calculator, local_feature_classifier,
                      patience, max_offspring_attempts, cut_and_splice_p,
                      selection_pressure, mutation_swap_fraction,
                      improvement_tol, seed):
    """Rust equivalent of ``run_single_particle_ga`` for the ACT case.

    Every member of the population is the same lattice with a different
    ordering, so only the symbol vectors cross the boundary: the geometry, the
    neighbour list and the masses are sent once. That is what removes the
    ``deepcopy`` per mutation and the neighbour-list rebuild per cut-and-splice
    that dominate the Python GA.

    Returns the same ``[best_energies, best_particle, energy_evaluations,
    population, stats]`` the Python returns.
    """
    from ase.data import atomic_masses

    reference = start_population[0]
    offsets, neighbours, _codes, elements = pack(reference)
    positions = np.asarray(reference.get_ase_atoms().get_positions(),
                           dtype=float).reshape(-1).tolist()
    masses = [float(atomic_masses[atomic_numbers[element]]) for element in elements]
    orderings = [pack(member)[2] for member in start_population]

    (final_orderings, final_energies, history, evaluations, generations,
     accepted, duplicates) = _rs.act_genetic_algorithm(
        offsets, neighbours, positions, masses, orderings,
        _coefficients(energy_calculator),
        [float(e) for e in environment_energies],
        int(patience), int(max_offspring_attempts), float(cut_and_splice_p),
        float(selection_pressure), float(mutation_swap_fraction[0]),
        float(mutation_swap_fraction[1]), float(improvement_tol), int(seed))

    energy_key = energy_calculator.get_energy_key()
    population = []
    for codes, energy in zip(final_orderings, final_energies):
        member = copy.deepcopy(reference)
        unpack(member, codes, elements)
        _restore(member, local_env_calculator, local_feature_classifier)
        member.set_energy(energy_key, energy)
        population.append(member)

    stats = dict(generations=generations, accepted=accepted,
                 duplicates=duplicates,
                 acceptance=accepted / max(generations, 1))
    return [history, population[0], evaluations, population, stats]


def local_optimization(start_particle, energy_calculator, environment_energies,
                       local_env_calculator, local_feature_classifier):
    """Rust equivalent of ``local_optimization`` for the ACT case."""
    offsets, neighbours, codes, elements = pack(start_particle)
    final_codes, accepted_energies, error = _rs.act_local_optimization(
        offsets, neighbours, codes, _coefficients(energy_calculator),
        [float(e) for e in environment_energies])

    particle = unpack(start_particle, final_codes, elements)
    _restore(particle, local_env_calculator, local_feature_classifier)
    energy_calculator.compute_energy(particle)
    return [particle, accepted_energies], error
