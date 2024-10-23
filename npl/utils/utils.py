from npl.core.nanoparticle import Nanoparticle
from npl.calculators import EMTCalculator
from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import learning_curve
from sklearn.metrics import mean_absolute_error,mean_squared_error
from sklearn.linear_model import BayesianRidge,Ridge
from sklearn.model_selection import ShuffleSplit
import numpy as np

def plot_learning_curves(X, y, n_atoms, estimator, n_splits=10, train_sizes=range(1, 401, 10),y_lim=(0,)):
    """
    Plots learning curves for a given estimator using Mean Absolute Error (MAE) as the scoring metric.

    Parameters:
    X (array-like): Feature matrix.
    y (array-like): Target vector.
    n_atoms (int): Number of atoms, used for normalizing the error.
    estimator (object): The estimator object implementing 'fit' and 'predict' methods.
    n_splits (int, optional): Number of re-shuffling & splitting iterations for cross-validation. Default is 10.
    train_sizes (iterable, optional): Relative or absolute numbers of training examples that will be used to generate the learning curve. Default is range(1, 401, 10).

    The function performs cross-validation to compute training and test scores, calculates the quartiles for the scores,
    and plots the learning curves with shaded areas representing the interquartile ranges.
    """
    cv = ShuffleSplit(n_splits=n_splits, train_size=train_sizes[-1], test_size=len(X) - train_sizes[-1], random_state=0)
    train_sizes, train_scores, test_scores = learning_curve(estimator, X=X, y=y, cv=cv, n_jobs=-1, train_sizes=train_sizes, scoring='neg_mean_absolute_error')

    train_scores = [np.quantile(train_scores, quartile, axis=1) / n_atoms * 1000 for quartile in [0.25, 0.50, 0.75]]
    test_scores = [np.quantile(test_scores, quartile, axis=1) / n_atoms * 1000 for quartile in [0.25, 0.50, 0.75]]

    train_q25, train_q50, train_q75 = train_scores
    test_q25, test_q50, test_q75 = test_scores

    plt.fill_between(train_sizes, -train_q25, -train_q75, alpha=0.3, label='Train IQR')
    plt.fill_between(train_sizes, -test_q25, -test_q75, alpha=0.3, label='Test IQR')
    plt.plot(train_sizes, -train_q50, '--', label='Train Median')
    plt.plot(train_sizes, -test_q50, '-', label='Test Median')
    plt.ylim(y_lim)
    plt.ylabel('MAE [meV / atom]')
    plt.xlabel('Training Set Size')
    plt.title('Learning Curves')
    plt.legend()
    plt.grid(True)
    plt.show()


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

def plot_cummulative_success_rate(energies: list, steps: list, figname: str = None):
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
    
    plt.ylim(0,100)
    plt.xlabel('Steps')
    plt.ylabel('Success Rate (%)')
    plt.title('Cumulative Success Rate')
    plt.grid(True)
    plt.show()
    
    if figname:
        plt.savefig(figname, dpi=200)