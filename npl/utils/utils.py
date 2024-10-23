from npl.core.nanoparticle import Nanoparticle
from npl.calculators import EMTCalculator
from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt

def get_reference_structure(particle: Nanoparticle, ase = None) -> Nanoparticle:
    if type(particle) is not Nanoparticle:
        p = Nanoparticle()
        p.add_atoms(particle)
        particle = p
    old_symbols = deepcopy(particle.atoms.atoms.symbols)
    new_symbols = ['Pt' for _ in particle.get_indices()]
    particle.transform_atoms(particle.get_indices(), new_symbols)
    EMTCalculator(relax_atoms=True).compute_energy(particle)
    particle.construct_neighbor_list()
    particle.transform_atoms(particle.get_indices(), old_symbols)
    if ase:
        return particle.get_ase_atoms()
    return particle

def plot_cummulative_success_rate(energies: list, steps: list, figname: str):
    """
    Plots the cumulative success rate based on given energies and steps, and saves the plot to a file.
    Parameters:
    energies (list): A list of energy values.
    steps (list): A list of step values corresponding to the energies.
    figname (str): The filename to save the plot.
    The function sorts the energies and steps, calculates the cumulative success rate, and plots it.
    The plot is saved as an image file with the specified filename.
    """
    energies, steps = zip(*sorted(zip(energies, steps)))
    min_energy = min(energies)
    max_steps = max(steps)
    success_rate = np.zeros(max_steps)
    
    percent = 0
    index = 0
    for step, energy in zip(steps, energies):
        if energy == min_energy:
            success_rate[step:] += 100 / len(energies)
            
    plt.plot(range(len(success_rate)), success_rate)
    plt.xlabel('Steps')
    plt.ylabel('Success Rate (%)')
    plt.title('Cumulative Success Rate')
    plt.grid(True)
    plt.savefig(figname, dpi=200)
    plt.show()