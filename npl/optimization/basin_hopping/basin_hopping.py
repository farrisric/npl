from npl.optimization.local_optimization import setup_local_optimization
from npl.optimization.local_optimization import update_atomic_features
import logging
import copy

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def run_basin_hopping(start_particle, energy_calculator, environment_energies, n_hopping_attempts,
                      n_hops,
                      local_feature_classifier=None):

    logging.info("Starting Basin Hopping simulation")
    logging.info("Hopping attempts: {}".format(n_hopping_attempts))
    logging.info("Hops: {}".format(n_hops))

    energy_key, local_env_calculator, local_feature_classifier, exchange_operator = \
        setup_local_optimization(start_particle, energy_calculator, environment_energies,
                                 local_feature_classifier)

    logging.info("Starting energy: {}".format(start_particle.get_energy(
        energy_calculator.get_energy_key())))

    start_energy = start_particle.get_energy(energy_key)
    lowest_energies = [(start_energy, 0)]
    best_particle = copy.deepcopy(start_particle)
    lowest_energy = start_energy

    flip_energy_list = []

    step = 0
    for i in range(n_hopping_attempts):
        while True:
            step += 1
            index1, index2 = exchange_operator.guided_exchange(start_particle)

            flip = exchange_operator.symbol1_exchange_energies[index1]
            + exchange_operator.symbol2_exchange_energies[index2]

            exchanged_indices = [index1, index2]

            start_particle, neighborhood = update_atomic_features(index1, index2,
                                                                  local_env_calculator,
                                                                  local_feature_classifier,
                                                                  start_particle)

            exchange_operator.update(start_particle, neighborhood, exchanged_indices)

            energy_calculator.compute_energy(start_particle)
            new_energy = start_particle.get_energy(energy_key)

            flip_energy_list.append((flip, new_energy-start_energy, exchanged_indices))

            if new_energy < start_energy:
                start_energy = new_energy
                if new_energy < lowest_energy:
                    lowest_energy = new_energy
                    lowest_energies.append((lowest_energy, step))

            else:
                if i % n_hopping_attempts == 0:
                    logging.info('Energy after local_opt: {:.3f}, lowest {:.3f}'.format(
                        start_energy,
                        lowest_energy))

                if lowest_energy == start_energy:
                    start_particle.swap_symbols([(index1, index2)])
                    best_particle = copy.deepcopy(start_particle)
                    best_particle.set_energy(energy_key, start_energy)

                    start_particle.swap_symbols([(index1, index2)])
                break

        for hop in range(n_hops):
            step += 1

            index1, index2 = exchange_operator.basin_hop_step(start_particle)

            exchanged_indices = [index1, index2]

            start_particle, neighborhood = update_atomic_features(index1, index2,
                                                                  local_env_calculator,
                                                                  local_feature_classifier,
                                                                  start_particle)
            exchange_operator.update(start_particle, neighborhood, exchanged_indices)

            energy_calculator.compute_energy(start_particle)
            new_energy = start_particle.get_energy(energy_key)

            start_energy = new_energy
    logging.info('Lowest energy: {:.3f}'.format(lowest_energy))
    lowest_energies.append((lowest_energy, step))

    return [best_particle, lowest_energies, flip_energy_list]
