import numpy as np
from typing import List
from itertools import combinations_with_replacement, product

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

    coefficeints = []
    feature_index_values = {Z : {} for Z in atomic_numbers}
    
    for Z_i in atomic_numbers:
        i = 0
        off_set_cn = len(bond_types) + n_species
        off_set_sublayer = len(bond_types)
        for cn_number in coordination_numbers:
            bond_energy = 0
            for bond_combination in compute_bond_combinations(n_species, cn_number):
                for n_bond, Z_j in zip(bond_combination, atomic_numbers):
                    
                    bond_type = tuple(sorted([Z_i, Z_j]))
                    bond_energy += n_bond * bond_types[bond_type]

                bond_energy += global_topological_coefficients[off_set_cn + cn_number]
                coefficeints.append(bond_energy)
                feature_index_values[Z_i][bond_combination] = i 
                i += 1
                if cn_number == 12:
                    bond_energy += global_topological_coefficients[off_set_sublayer]
                    coefficeints.append(bond_energy)
                    bond_combination = [x for x in bond_combination] + [1]
                    feature_index_values[Z_i][tuple(bond_combination)] = i 
                    i += 1
            
        off_set_sublayer += 1
        off_set_cn += 13

    return coefficeints, feature_index_values

def compute_bond_combinations(n_bond_types, cn_number):
    n_bonds = [x for x in reversed(range(13))]

    combinations = []
    for pair in product(n_bonds, repeat=n_bond_types):
        if sum(pair) == cn_number:
            combinations.append(pair)

    return combinations

    
