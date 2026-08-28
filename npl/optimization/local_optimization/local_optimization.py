from npl.descriptors.local_environment_calculator import NeighborCountingEnvironmentCalculator
from npl.optimization.local_optimization.guided_exchange_operator import GuidedExchangeOperator
from npl.optimization.local_optimization.act_exchange_operator import ACTExchangeOperator
from npl.descriptors.local_environment_feature_classifier import TopologicalEnvironmentClassifier
from npl.descriptors.global_feature_classifier import ExtendedTopologicalFeaturesClassifier
from npl.optimization.local_optimization.etop_exchange_operator import ETOPExchangeOperator
import copy

EXCHANGE_OPERATORS = {'TOP': GuidedExchangeOperator, 'ACT': ACTExchangeOperator,
                      'ETOP': ETOPExchangeOperator}


def setup_local_optimization(start_particle, energy_calculator, environment_energies,
                             local_feature_classifier=None, model='TOP'):
    symbols = start_particle.get_all_symbols()
    if model == 'ETOP':
        if not hasattr(energy_calculator, 'get_coefficients'):
            raise ValueError(
                "ETOP guided exchange requires a linear energy model exposing "
                "get_coefficients() (e.g. BayesianRRCalculator or TOPCalculator).")
        if local_feature_classifier is None:
            local_feature_classifier = ExtendedTopologicalFeaturesClassifier(symbols)
        local_feature_classifier.compute_feature_vector(start_particle)
        energy_calculator.compute_energy(start_particle)
        exchange_operator = ETOPExchangeOperator(
            energy_calculator.get_coefficients(), local_feature_classifier)
        exchange_operator.bind_particle(start_particle)
        energy_key = energy_calculator.get_energy_key()
        return energy_key, None, local_feature_classifier, exchange_operator
    local_env_calculator = NeighborCountingEnvironmentCalculator(symbols)
    if local_feature_classifier is None:
        local_feature_classifier = TopologicalEnvironmentClassifier(local_env_calculator, symbols)

    local_env_calculator.compute_local_environments(start_particle)
    local_feature_classifier.compute_feature_vector(start_particle)

    feature_key = local_feature_classifier.get_feature_key()
    energy_calculator.compute_energy(start_particle)

    try:
        exchange_operator_cls = EXCHANGE_OPERATORS[model]
    except KeyError:
        raise ValueError(
            f"Unknown model {model!r}; choose from {sorted(EXCHANGE_OPERATORS)}")
    exchange_operator = exchange_operator_cls(environment_energies, feature_key)
    exchange_operator.bind_particle(start_particle)

    energy_key = energy_calculator.get_energy_key()

    return energy_key, local_env_calculator, local_feature_classifier, exchange_operator


def update_atomic_features(index1, index2, local_env_calculator, local_feature_classifier,
                           particle):
    neighborhood = {index1, index2}
    neighborhood = neighborhood.union(particle.neighbor_list[index1])
    neighborhood = neighborhood.union(particle.neighbor_list[index2])

    for x in neighborhood:
        neighborhood = neighborhood.union(particle.neighbor_list[x])

    for index in neighborhood:
        local_env_calculator.compute_local_environment(particle, index)
        local_feature_classifier.compute_atom_feature(particle, index)

    local_feature_classifier.compute_feature_vector(particle, recompute_atom_features=False)

    return particle, neighborhood


def update_atomic_features_etop(index1, index2, feature_classifier, particle):
    neighborhood = {index1, index2}
    neighborhood = neighborhood.union(particle.neighbor_list[index1])
    neighborhood = neighborhood.union(particle.neighbor_list[index2])

    old_atom_features, change = feature_classifier.update_feature_vector(particle, neighborhood)

    return particle, neighborhood, old_atom_features, change


def local_optimization(start_particle, energy_calculator, environment_energies,
                       local_feature_classifier=None, model='TOP',
                       force_python=False):
    """Guided-exchange descent to a local minimum of the chemical ordering.

    ``model='ACT'`` on a two-element particle with a linear TEC energy model
    runs inside the ``npl_act_rs`` extension when it is installed. The descent
    is deterministic, so the two paths agree exactly except where equal exchange
    energies leave the choice of atom open. ``force_python=True`` keeps the pure
    Python path.
    """
    from npl.optimization import act_fast  # lazy: npl.optimization imports this module
    from npl.optimization.basin_hopping.basin_hopping import _act_descriptors

    if not force_python and act_fast.usable(start_particle, energy_calculator,
                                            environment_energies, model):
        env_calculator, feature_classifier = _act_descriptors(
            start_particle, local_feature_classifier)
        env_calculator.compute_local_environments(start_particle)
        feature_classifier.compute_feature_vector(start_particle)
        result, _error = act_fast.local_optimization(
            copy.deepcopy(start_particle), energy_calculator, environment_energies,
            env_calculator, feature_classifier)
        return result

    energy_key, local_env_calculator, local_feature_classifier, exchange_operator = \
        setup_local_optimization(start_particle, energy_calculator,
                                 environment_energies, local_feature_classifier, model=model)

    start_energy = start_particle.get_energy(energy_key)
    accepted_energies = [(start_energy, 0)]

    step = 0
    while True:
        step += 2
        index1, index2, _flip = exchange_operator.guided_exchange(start_particle)
        exchanged_indices = [index1, index2]

        if local_env_calculator is None:
            start_particle, neighborhood, old_atom_features, change = \
                update_atomic_features_etop(index1, index2, local_feature_classifier,
                                            start_particle)
        else:
            start_particle, neighborhood = update_atomic_features(
                index1, index2, local_env_calculator, local_feature_classifier, start_particle)
        exchange_operator.update(start_particle, neighborhood, exchanged_indices)

        energy_calculator.compute_energy(start_particle)
        new_energy = start_particle.get_energy(energy_key)

        if new_energy < start_energy:
            start_energy = new_energy
            accepted_energies.append((new_energy, step))
        else:
            accepted_energies.append((start_energy, step))

            # roll back last exchange and make sure features and environments are up-to-date
            start_particle.swap_symbols([(index1, index2)])
            start_particle.set_energy(energy_key, start_energy)
            if local_env_calculator is None:
                local_feature_classifier.downgrade_feature_vector(
                    start_particle, neighborhood, old_atom_features, change)
            else:
                for index in neighborhood:
                    local_env_calculator.compute_local_environment(start_particle, index)
                    local_feature_classifier.compute_atom_feature(start_particle, index)

            break

    return [start_particle, accepted_energies]
