from itertools import combinations_with_replacement, product
from typing import List

import numpy as np


def get_combinations(element_list : List, combination_size : int) -> List[List]:
    "Returns a list of combinations of a list with given size"
    return [x for x in combinations_with_replacement(element_list, combination_size)]

def compute_bond_types(coefficients, atomic_numbers):
    bond_type_and_coefficient = {}
    for bond_type, coefficient in zip(combinations_with_replacement(atomic_numbers, 2), coefficients):
        bond_type_and_coefficient[bond_type] = coefficient
    return bond_type_and_coefficient

def compute_coefficients_for_linear_topological_model(global_topological_coefficients, atomic_numbers):
    atomic_numbers = sorted(atomic_numbers)
    n_species = len(atomic_numbers)
    coordination_numbers = list(range(13))
    bond_types = compute_bond_types(global_topological_coefficients, atomic_numbers)

    coefficients = []
    total_energies = []
    feature_index_values = {Z : {} for Z in atomic_numbers}

    i = 0
    for Z_i in atomic_numbers:
        off_set_cn = len(bond_types) + n_species
        off_set_sublayer = len(bond_types)
        for cn_number in coordination_numbers:
            for bond_combination in compute_bond_combinations(n_species, cn_number):
                env_energy = 0
                total_energy = 0
                for n_bond, Z_j in zip(bond_combination, atomic_numbers):
                    bond_type = tuple(sorted([Z_i, Z_j]))

                    env_energy += ((n_bond/2) * bond_types[bond_type])
                    total_energy += (n_bond * bond_types[bond_type])

                env_energy += global_topological_coefficients[off_set_cn + cn_number]
                total_energy += global_topological_coefficients[off_set_cn + cn_number]

                coefficients.append(env_energy)
                total_energies.append(total_energy)
                feature_index_values[Z_i][bond_combination] = i 
                i += 1
                # if cn_number == 12:
                #     env_energy += (global_topological_coefficients[off_set_sublayer])
                #     total_energies += (global_topological_coefficients[off_set_sublayer])
                #     coefficients.append(env_energy)
                #     total_energies.append(total_energies)
                #     bond_combination = [x for x in bond_combination] + [1]
                #     feature_index_values[Z_i][tuple(bond_combination)] = i 
                #     i += 1
            
        off_set_sublayer += 1
        off_set_cn += 13

    return coefficients, total_energies, feature_index_values

def compute_bond_combinations(n_bond_types, cn_number):
    n_bonds = [x for x in reversed(range(13))]

    combinations = []
    for pair in product(n_bonds, repeat=n_bond_types):
        if sum(pair) == cn_number:
            combinations.append(pair)

    return combinations

    
