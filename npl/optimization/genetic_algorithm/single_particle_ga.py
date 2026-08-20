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
    if max_energy == min_energy:
        return UNIFORM
    normalized_energy = (particle.get_energy(energy_key) - min_energy) / (max_energy - min_energy)
    return np.exp(-3 * normalized_energy)


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
                          max_offspring_attempts=50, cut_and_splice_p=0.4):
    """Steady-state genetic algorithm over the orderings of one particle.

    ``model`` selects the exchange operator used to locally optimize every
    offspring. It defaults to 'ACT' because the 'TOP' GuidedExchangeOperator is
    fully deterministic - both its guided step and its kick - so every offspring
    descends into the same limit cycle and the population collapses.

    Convergence is counted with ``improvement_tol`` rather than exact float
    equality: an improvement of 1e-12 is float noise, and treating it as
    progress kept the loop alive indefinitely.
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

    while unsuccessful_gens < unsuccessful_gens_for_convergence:
        min_energy = cur_population[0].get_energy(energy_key)
        max_energy = cur_population[-1].get_energy(energy_key)
        generation += 1

        fitness_values = np.array(
            [compute_fitness(p, min_energy, max_energy, energy_key) for p in cur_population])
        fitness_values /= np.sum(fitness_values)

        new_offspring = None
        for _attempt in range(max_offspring_attempts):
            if np.random.rand() < cut_and_splice_p and len(cur_population) > 1:
                parent1, parent2 = np.random.choice(
                    cur_population, 2, replace=False, p=fitness_values)
                candidate = cut_and_splice_operator.cut_and_splice(parent1, parent2)
            else:
                parent = np.random.choice(cur_population, 1, p=fitness_values)[0]
                candidate = exchange_operator.random_exchange(parent)

            candidate, energies = local_optimization(
                candidate, energy_calculator, environment_energies,
                model=model)
            energy_evaluations += energies[-1][1]

            if ordering_key(candidate) not in seen:
                new_offspring = candidate
                break

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

        best_energy = cur_population[0].get_energy(energy_key)
        if best_energy < best_energies[-1][0] - improvement_tol:
            unsuccessful_gens = 0
            best_energies.append((best_energy, energy_evaluations))
            logger.info('generation %d: new best energy %.6f', generation, best_energy)
        else:
            unsuccessful_gens += 1

    return [best_energies, cur_population[0], energy_evaluations]
