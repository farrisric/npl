Examples
========

This section provides practical examples of how to use the Nanoparticle Library (NPL) to train models, extract features, and optimize nanoparticle structures.

TOP Energy Evaluation
---------------------

The flagship example shows how to evaluate the topological (TOP) energy of a bimetallic nanoparticle using pretrained coefficients shipped with NPL. See ``examples/top_energy_evaluation.py`` in the repository.

Training a Topological Descriptor Model
----------------------------------------

For a detailed tutorial on training a model using Topological Descriptors, please refer to the following document:

.. toctree::
   :maxdepth: 1

   examples/train_top.rst
   examples/multimet_go.rst

Basin hopping
-------------

``examples/basin_hopping_top.py`` and ``examples/basin_hopping_act.py`` optimize the chemical ordering of a two-metal particle with guided-exchange basin hopping (``model="TOP"`` and ``model="ACT"``).

.. note::

   Canonical Monte Carlo over chemical ordering is part of NPL (``npl.monte_carlo``, see the multimetallic example).
   Grand-canonical Monte Carlo with machine-learning potentials lives in the companion package mcpy (https://github.com/farrisric/mcpy).
