Loading a Pretrained TOP Model and Performing Chemical Ordering Optimization
==========================================================================

This tutorial demonstrates how to load a pretrained TOP model and perform a chemical ordering optimization using the `npl` library.

First, we need to import the necessary modules and initialize the `TOPCalculator` with the `ExtendedTopologicalFeaturesClassifier`.

.. code-block:: python

    from npl.descriptors import ExtendedTopologicalFeaturesClassifier
    from npl.calculators import TOPCalculator

    calc = TOPCalculator('ETOP', stoichiometry='Pt151Cu50',
                         feature_classifier=ExtendedTopologicalFeaturesClassifier)

    etop = calc.get_feature_classifier()

When the `TOPCalculator` is initialized, it will load the topological parameters for the given stoichiometry. You should see output similar to the following:

.. code-block:: text

    INFO - Loading top parameters of Pt151Cu50
    INFO - Parameters obtained from reference: L. Vega Mater. Adv., 2021, 2, 6589-6602
    INFO - Parameters loaded successfully
    INFO - Parameters: 
    {'CuPt': -25.0, 'Cu(cn=6)': 267.0, 'Cu(cn=7)': 342.0, 'Cu(cn=8)': 372.0, 'Cu(cn=9)': 372.0}

Next, we will run 20 Monte Carlo simulations to optimize the chemical ordering.

    .. code-block:: python

        from npl.monte_carlo import run_monte_carlo
        from npl.core import Nanoparticle

        beta = 250
        max_steps = 10000

        energy_calculator = calc
        feature_classifier = etop

        energies_MC, steps_MC = [], []
        for i in range(10):
            start_particle = Nanoparticle()
            start_particle.truncated_octahedron(7, 2, {'Pt': 151, 'Cu': 50})
            best_particle, accepted_energies = run_monte_carlo(beta, max_steps,
                                                               start_particle,
                                                               energy_calculator,
                                                               feature_classifier)
            min_energy, min_step = min(accepted_energies, key=lambda x: x[0])
            energies_MC.append(min_energy)
            steps_MC.append(min_step)
            if min_energy <= min(energies_MC):
                global_minimum = best_particle


Finally, we evaluate the results of our Monte Carlo simulations by looking at the cumulative success rate plot and visualizing the global minimum chemical ordering. We also plot its resulting topological descriptor.

.. code-block:: python

    import matplotlib.pyplot as plt
    from npl.visualization import plot_topological_descriptor

    # Plot cumulative success rate
    plt.figure()
    plt.plot(range(1, len(energies_MC) + 1), sorted(energies_MC))
    plt.xlabel('Simulation Run')
    plt.ylabel('Minimum Energy')
    plt.title('Cumulative Success Rate')
    plt.show()

    # Visualize global minimum chemical ordering
    global_minimum.visualize()

    # Plot topological descriptor of the global minimum
    plot_topological_descriptor(global_minimum)
