# NanoParticleLibrary (NPL)

NPL (NanoParticleLibrary) is a comprehensive wrapper around the popular ASE (Atomic Simulation Environment) library, designed to facilitate the manipulation and optimization of nanoparticles. It is particularly tailored for optimizing the chemical ordering in bimetallic nanoparticles. 

## Features

NPL offers a variety of features to support nanoparticle research and optimization:

- **Flexible Energy Pipeline**: Supports approximate energy models for efficient computations.
- **Local Optimization**: Optimizes the ordering of particles locally.
- **Pre-built Optimization Algorithms**: Includes algorithms such as Markov Chain Monte Carlo (MCMC) and Genetic Algorithms (GA).
- **Basin-Hopping**: Utilizes local optimization for both relaxation and perturbation steps.
- **Shape Optimization**: Provides built-in procedures for optimizing the shape of fixed lattice structures (e.g., FCC, icosahedral, decahedral).
- **Memory Efficient Storage**: Efficiently stores large training sets to save memory.

## Installation

To install NPL, ensure you have the following dependencies:

- Python 3.7+
- Atomic Simulation Environment (ASE) >= 3.21
- scikit-learn
- sortedcontainers (by Grant Jenks)
