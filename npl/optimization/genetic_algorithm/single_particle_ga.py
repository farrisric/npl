import logging

import numpy as np

from npl.optimization.genetic_algorithm.cut_and_splice_operator import CutAndSpliceOperator
from npl.optimization.genetic_algorithm.exchange_operator import ExchangeOperator
from npl.optimization.local_optimization.local_optimization import local_optimization

logger = logging.getLogger(__name__)

# A population whose members all share one energy carries no selection signal.
# Returning zero fitness for every member used to break the generation loop,
# which looked exactly like convergence; a uniform fitness keeps the search
# alive and lets the convergence counter decide when to stop.
UNIFORM = 1.0


def compute_fitness(particle, min_energy, max_energy, energy_key):
    """Legacy fitness: energy normalized by the population's own spread.

    Kept for callers that want it, but note the coupling it introduces - the
    same set of energy gaps means something different at 7 eV of spread than at
    0.1 eV, so selection pressure drifts with the population's diversity
    instead of tracking a physical scale. Prefer `rank_fitness`.
    """
    if max_energy == min_energy:
        return UNIFORM
    normalized_energy = (particle.get_energy(energy_key) - min_energy) / (max_energy - min_energy)
    return np.exp(-3 * normalized_energy)


def rank_fitness(n_members, pressure):
    """Selection weights from rank alone, for a population sorted best-first.

    Scale-free by construction: the best member always gets weight 1 and the
    worst exp(-pressure), whatever the energies happen to be. `pressure` is
    then a real knob - 0 is a random walk, large values are near-greedy.
    """
    if n_members == 1:
        return np.array([UNIFORM])
    ranks = np.arange(n_members) / (n_members - 1)
    return np.exp(-pressure * ranks)


def mutation_size(n_atoms, swap_fraction_range, rng=np.random):
    """How many swaps one mutation applies.

    A single swap is close to useless here: the parent has already been
    descended to a local minimum of the surrogate, so by definition no single
    swap improves it, and the offspring is discarded. Mutations therefore have
    to be real kicks, sized as a fraction of the particle.
    """
    low, high = swap_fraction_range
    lo = max(1, int(round(low * n_atoms)))
    hi = max(lo, int(round(high * n_atoms)))
    return int(rng.randint(lo, hi + 1))


def ordering_key(particle):
    """Identity of an ordering, for de-duplication.

    Energies are the wrong key: on a symmetric particle many distinct orderings
    are exactly degenerate, so an energy-based test throws away the diversity a
    genetic algorithm depends on.
    """
    return tuple(particle.get_symbols())


def run_single_particle_ga(start_population, unsuccessful_gens_for_convergence,
                          energy_calculator, local_env_calculator,
                          local_feature_classifier, environment_energies,
                          model='ACT', improvement_tol=1e-6,
                          max_offspring_attempts=50, cut_and_splice_p=0.4,
                          selection_pressure=3.0,
                          mutation_swap_fraction=(0.01, 0.06)):
    """Steady-state genetic algorithm over the orderings of one particle.

    ``model`` selects the exchange operator used to locally optimize every
    offspring. It defaults to 'ACT' because the 'TOP' GuidedExchangeOperator is
    fully deterministic - both its guided step and its kick - so every offspring
    descends into the same limit cycle and the population collapses.

    Convergence is counted with ``improvement_tol`` rather than exact float
    equality: an improvement of 1e-12 is float noise, and treating it as
    progress kept the loop alive indefinitely.

    Selection is rank-based with ``selection_pressure``, so it does not drift
    with the population's spread, and mutations swap
    ``mutation_swap_fraction`` of the atoms rather than the one atom a
    geometric draw used to give.

    Returns ``[best_energies, best_particle, energy_evaluations, population,
    stats]``; ``stats`` reports generations, how many offspring survived
    selection, and how many were rejected as duplicates, which is what tells
    you whether the settings are wasting the budget.
    """
    unsuccessful_gens = 0
    energy_key = energy_calculator.get_energy_key()

    for p in start_population:
        local_env_calculator.compute_local_environments(p)
        local_feature_classifier.compute_feature_vector(p)
        energy_calculator.compute_energy(p)

    exchange_operator = ExchangeOperator(0.5)
    cut_and_splice_operator = CutAndSpliceOperator()

    cur_population = list(start_population)
    generation = 0
    best_energies = []
    energy_evaluations = len(cur_population)

    cur_population.sort(key=lambda x: x.get_energy(energy_key))
    best_energies.append((cur_population[0].get_energy(energy_key), energy_evaluations))
    seen = {ordering_key(p) for p in cur_population}
    n_accepted, n_duplicates = 0, 0
    n_atoms = len(cur_population[0].get_symbols())

    while unsuccessful_gens < unsuccessful_gens_for_convergence:
        generation += 1

        # population is kept sorted best-first, so rank weights need no energies
        fitness_values = rank_fitness(len(cur_population), selection_pressure)
        fitness_values /= np.sum(fitness_values)

        new_offspring = None
        for _attempt in range(max_offspring_attempts):
            if np.random.rand() < cut_and_splice_p and len(cur_population) > 1:
                parent1, parent2 = np.random.choice(
                    cur_population, 2, replace=False, p=fitness_values)
                candidate = cut_and_splice_operator.cut_and_splice(parent1, parent2)
            else:
                parent = np.random.choice(cur_population, 1, p=fitness_values)[0]
                candidate = exchange_operator.random_exchange(
                    parent,
                    n_exchanges=mutation_size(n_atoms, mutation_swap_fraction))

            candidate, energies = local_optimization(
                candidate, energy_calculator, environment_energies,
                model=model)
            energy_evaluations += energies[-1][1]

            if ordering_key(candidate) not in seen:
                new_offspring = candidate
                break
            n_duplicates += 1

        if new_offspring is None:
            # every attempt reproduced an ordering already in the population:
            # the population has collapsed and more generations cannot help
            logger.info('GA stopped at generation %d: no new ordering in %d '
                        'attempts', generation, max_offspring_attempts)
            break

        seen.add(ordering_key(new_offspring))
        cur_population.append(new_offspring)
        cur_population.sort(key=lambda x: x.get_energy(energy_key))
        dropped = cur_population.pop()
        seen.discard(ordering_key(dropped))
        if dropped is not new_offspring:
            n_accepted += 1

        best_energy = cur_population[0].get_energy(energy_key)
        if best_energy < best_energies[-1][0] - improvement_tol:
            unsuccessful_gens = 0
            best_energies.append((best_energy, energy_evaluations))
            logger.info('generation %d: new best energy %.6f', generation, best_energy)
        else:
            unsuccessful_gens += 1

    stats = dict(generations=generation, accepted=n_accepted,
                 duplicates=n_duplicates,
                 acceptance=n_accepted / max(generation, 1))
    logger.info('GA finished: %d generations, %d offspring accepted (%.0f%%), '
                '%d duplicates rejected', generation, n_accepted,
                100 * stats['acceptance'], n_duplicates)
    return [best_energies, cur_population[0], energy_evaluations,
            cur_population, stats]
