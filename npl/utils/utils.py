import numpy as np
from itertools import combinations_with_replacement, product

def compute_bond_types(coefficients, symbols):

    bond_type_and_coefficient = {}
    for bond_type, coefficient in zip(combinations_with_replacement(symbols, 2), coefficients):
        bond_type_and_coefficient[bond_type] = coefficient

    return bond_type_and_coefficient

def compute_coefficients_for_linear_topological_model(global_topological_coefficients, symbols):
    symbols = sorted(symbols)
    n_symbol = len(symbols)
    coordination_numbers = list(range(13))
    bond_types = compute_bond_types(global_topological_coefficients, symbols)

    coefficeints = []
    feature_index_values = {symbol : {} for symbol in symbols}
    
    for symbol_a in symbols:
        i = 0
        off_set_cn = len(bond_types) + n_symbol
        off_set_sublayer = len(bond_types)
        for cn_number in coordination_numbers:
            bond_energy = 0
            for bond_combination in compute_bond_combinations(n_symbol, cn_number):
                for n_bond, symbol_b in zip(bond_combination, symbols):
                    
                    bond_type = tuple(sorted([symbol_a, symbol_b]))
                    bond_energy += n_bond * bond_types[bond_type]

                bond_energy += global_topological_coefficients[off_set_cn + cn_number]
                coefficeints.append(bond_energy)
                feature_index_values[symbol_a][bond_combination] = i 
                i += 1
                if cn_number == 12:
                    bond_energy += global_topological_coefficients[off_set_sublayer]
                    coefficeints.append(bond_energy)
                    bond_combination = [x for x in bond_combination] + [1]
                    feature_index_values[symbol_a][tuple(bond_combination)] = i 
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



if __name__ == "__main__":
    import random 

    global_top = [random.random() for _ in range(48)]
    symbols = ['Au','Pt','Ni']
    # a = compute_bond_types(global_top, symbols)
    coefficeints, feature_index_values = compute_coefficients_for_linear_topological_model(global_top, symbols)
    print(feature_index_values['Pt'][(0, 1, 11, 1)])
    
