import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import pandas as pd
from sklearn.model_selection import learning_curve, ShuffleSplit
from ase.visualize.plot import plot_atoms


def plot_cummulative_success_rate(energies: list, steps: list, figname: str = None):
    """
    Plots the cumulative success rate based on given energies and steps, and saves the plot to a
    file.
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
    for step, energy in zip(steps, energies):
        if energy == min_energy:
            success_rate[step:] += 100 / len(energies)

    plt.plot(range(len(success_rate)), success_rate)

    plt.ylim(0, 100)
    plt.xlabel('Steps')
    plt.ylabel('Success Rate (%)')
    plt.title('Cumulative Success Rate')
    plt.grid(True)
    plt.show()

    if figname:
        plt.savefig(figname, dpi=200)


def get_surface_core_indices(particle):
    surface_indices = particle.get_atom_indices_from_coordination_number(range(12))
    core_indices = particle.get_atom_indices_from_coordination_number([12])
    return surface_indices, core_indices


def separate_layers(particle, upper_layer, indices):
    new_upper = []
    new_indices = []
    for idx in indices:
        neighbor_list = particle.neighbor_list[idx]
        cacca = False
        for x in neighbor_list:
            if x in upper_layer:
                new_upper.append(idx)
                cacca = True
                break
        if not cacca:
            new_indices.append(idx)
    return new_upper, new_indices


def get_layers_indices(particle):
    surface_indices, core_indices = get_surface_core_indices(particle)
    upper_layer = surface_indices
    indices = core_indices
    layers = [surface_indices]
    while indices:
        upper_layer, indices = separate_layers(particle, upper_layer, indices)
        layers.append(upper_layer)
    return layers


def get_concentration_per_layer(particle, symbols):
    layers = get_layers_indices(particle)
    layers_n_atoms = [len(x) for x in layers]
    concentration_layers = []
    for layer in layers:
        layer_symbols = particle.get_symbols(layer)
        layer_concentration = [layer_symbols.count(symbol) / len(layer) for symbol in symbols]
        concentration_layers.append(layer_concentration)
    return concentration_layers, layers_n_atoms


def plot_elemental_concentration_per_layer(particle):
    from ase.data import colors
    from ase.data import atomic_numbers
    """
    Plots the elemental concentration per layer for a given particle.
    Parameters:
    particle (object): The particle object containing the data.
    symbols (list of str): List of element symbols to be plotted.
    cmap (str or Colormap): Colormap to be used for the plot.
    colors (list of str or list of tuple): List of colors to be used for the plot.
    Returns:
    None
    """
    symbols = particle.get_all_symbols()
    numbers = [atomic_numbers[symbol] for symbol in symbols]
    atom_colors = [colors.jmol_colors[number] for number in numbers]
    concentration_layers, layers_n_atoms = get_concentration_per_layer(particle, symbols)
    layer_conc = pd.DataFrame(concentration_layers, columns=symbols)
    cmap = LinearSegmentedColormap.from_list("", atom_colors)

    fig, ax = plt.subplots(figsize=(12, 8))
    layer_conc.plot(kind='bar',
                    stacked=True,
                    colormap=cmap,
                    width=1,
                    edgecolor='k',
                    ax=ax)

    for i, (concentration, n_atoms) in enumerate(zip(concentration_layers, layers_n_atoms)):
        for j, (symbol, conc) in enumerate(zip(symbols, concentration)):
            if conc < 0.01:
                continue
            ax.text(i, sum(concentration[:j]) + conc / 2, f'{conc * 100:.1f}%',
                    ha='center', va='center', color='black', fontsize=10, fontweight='bold')

    ax.set_xlabel('Shell Number (0=Surface, {}=Core)'.format(len(concentration_layers) - 1),
                  fontsize=14)
    ax.set_ylabel('Shell Composition', fontsize=14)
    ax.set_xlim(-0.5, len(layer_conc) - 0.5)
    ax.set_ylim(0, 1)
    ax.set_title('Elemental Concentration per Layer', fontsize=16)
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)
    ax.legend(title='Elements', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=12)
    plt.tight_layout()
    plt.show()


def plot_learning_curves(X, y, n_atoms, estimator, n_splits=10,
                         train_sizes=range(1, 401, 10), y_lim=None,
                         filename=None):
    """
    Plots learning curves for a given estimator using Mean Absolute Error (MAE)
    as the scoring metric.

    Parameters:
    X (array-like): Feature matrix.
    y (array-like): Target vector.
    n_atoms (int): Number of atoms, used for normalizing the error.
    estimator (object): The estimator object implementing 'fit' and 'predict' methods.
    n_splits (int, optional): Number of re-shuffling & splitting iterations for cross-validation.
    Default is 10.
    train_sizes (iterable, optional): Numbers of training examples used to generate the learning
    curve. Default is range(1, 401, 10).

    The function performs cross-validation to compute training and test scores, calculates the
    quartiles for the scores, and plots the learning curves with shaded areas representing the
    interquartile ranges.
    """

    # Cross-validation setup
    cv = ShuffleSplit(n_splits=n_splits, train_size=train_sizes[-1],
                      test_size=len(X) - train_sizes[-1])

    # Calculate learning curves
    train_sizes, train_scores, test_scores = learning_curve(
        estimator, X=X, y=y, cv=cv, n_jobs=-1,
        train_sizes=train_sizes, scoring='neg_mean_absolute_error'
    )

    # Calculate quartiles for training and test scores
    train_scores = [np.quantile(train_scores, quartile, axis=1) for quartile in [0.25, 0.50, 0.75]]
    test_scores = [np.quantile(test_scores, quartile, axis=1) for quartile in [0.25, 0.50, 0.75]]

    train_q25, train_q50, train_q75 = train_scores
    test_q25, test_q50, test_q75 = test_scores

    # Plotting the learning curves
    plt.figure(figsize=(10, 6))
    plt.fill_between(train_sizes, -train_q25, -train_q75,
                     alpha=0.3, label='Train IQR', color='lightblue')
    plt.fill_between(train_sizes, -test_q25, -test_q75,
                     alpha=0.3, label='Test IQR', color='lightgreen')

    plt.plot(train_sizes, -train_q50, '--', label='Train Median', color='blue')
    plt.plot(train_sizes, -test_q50, '-', label='Test Median', color='green')
    if y_lim:
        plt.ylim(y_lim)
    plt.ylabel('MAE [meV / atom]', fontsize=12)
    plt.xlabel('Training Set Size', fontsize=12)
    plt.title('Learning Curves', fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(fontsize=10)

    # Save the figure if a filename is provided
    if filename:
        plt.savefig(filename, dpi=200)


def plot_parted_particle(atoms, separation=3):
    atoms = atoms.get_ase_atoms()
    atoms.center()
    atoms1 = atoms[[a.index for a in atoms if a.position[2] < atoms.get_cell()[2][2]/2 + 1.0]]
    atoms2 = atoms[[a.index for a in atoms if a.position[2] > atoms.get_cell()[2][2]/2 + 1.0]]
    atoms1.translate((0., 0., -separation))
    atoms2.translate((0., 0., separation))
    atoms = atoms1 + atoms2
    plot_atoms(atoms, rotation=('0x,75y,0z'))
    plt.axis('off')
