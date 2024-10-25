Examples
========

Basic Usage Example
-------------------

This example demonstrates setting up a bimetallic nanoparticle and optimizing it with Monte Carlo:

.. code-block:: python

   from nplib import Nanoparticle
   from nplib.optimization import MonteCarlo
   
   # Initialize a bimetallic nanoparticle
   nanoparticle = Nanoparticle(elements=["Ag", "Pd"], size=50)

   # Run Monte Carlo optimization
   optimizer = MonteCarlo(nanoparticle)
   optimized_structure = optimizer.run()

Training a Surrogate Model
--------------------------

This example shows how to train a surrogate model with topological descriptors:

.. code-block:: python

   from nplib.models import SurrogateModel
   from nplib.descriptors import TopologicalDescriptor

   # Initialize descriptors and surrogate model
   descriptors = TopologicalDescriptor(nanoparticle)
   model = SurrogateModel(descriptors)

   # Train the model
   model.train()
