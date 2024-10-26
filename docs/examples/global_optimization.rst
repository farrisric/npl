Multimetal Optimization Tutorial
================================

This tutorial demonstrates how to use the NPL (NanoParticleLibrary) for optimizing the chemical ordering in multimetallic nanoparticles.

.. contents:: Table of Contents
   :depth: 2
   :local:

Prerequisites
-------------
- Basic understanding of optimization techniques.
- Familiarity with Python programming.
- Knowledge of numerical libraries such as NumPy.

Getting Started
---------------
To get started, ensure you have the NPL library installed. You can install it using pip:

.. code-block:: bash

    pip install npl

Importing Necessary Modules
---------------------------
First, import the necessary modules from the NPL library:

.. code-block:: python

    import npl.core.nanoparticle as NP
    import npl.optimization.go_search as GOS
    from npl.optimization.local_optimization.local_optimization import local_optimization

Creating a Nanoparticle
-----------------------
Create a nanoparticle object using the `npl.core.nanoparticle` module:

.. code-block:: python

    # Example code to create a nanoparticle
    nanoparticle = NP.Nanoparticle('path/to/structure/file')

Setting Up the Optimization
---------------------------
Set up the optimization parameters and initialize the optimization process:

.. code-block:: python

    # Example code to set up optimization
    optimizer = GOS.GlobalOptimizationSearch(nanoparticle)
    optimizer.setup_parameters(param1=value1, param2=value2)

Running the Optimization
------------------------
Run the optimization process and monitor the progress:

.. code-block:: python

    # Example code to run optimization
    results = optimizer.run()
    print(results)

Analyzing the Results
---------------------
Analyze the optimization results and visualize the optimized structure:

.. code-block:: python

    # Example code to analyze results
    optimized_structure = results.get_optimized_structure()
    optimized_structure.visualize()

Conclusion
----------
This tutorial provided a step-by-step guide on how to use the NPL library for optimizing the chemical ordering in multimetallic nanoparticles. For more detailed examples and advanced usage, refer to the official documentation.

References
----------
If you use this code, please cite our papers:

.. code-block:: bibtex

    @article{10.1063/5.0193848,
    author = {Farris, Riccardo and Merinov, Boris V. and Bruix, Albert and Neyman, Konstantin M.},
    title = "{Effects of Zr dopants on properties of PtNi nanoparticles for ORR catalysis: A DFT modeling}",
    journal = {The Journal of Chemical Physics},
    volume = {160},
    number = {12},
    pages = {124706},
    year = {2024},
    issn = {0021-9606},
    doi = {10.1063/5.0193848},
    url = {https://doi.org/10.1063/5.0193848},
}