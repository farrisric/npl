<!-- [![GitHub release](https://img.shields.io/github/release/yourusername/npl.svg)](https://GitHub.com/yourusername/npl/releases/) -->

[![License](https://img.shields.io/badge/License-MIT%202.0-blue.svg)](https://opensource.org/licenses/mit)
[![GitHub issues](https://img.shields.io/github/issues/farrisric/nplib.svg)](https://GitHub.com/farrisric/NPlib/issues)
[![Documentation Status](https://readthedocs.org/projects/nplib/badge/)](https://nplib.readthedocs.io/en/latest/index.html)

# <span style="font-size:larger;">NanoParticleLibrary (NPL)</span>

![NPL Logo](https://github.com/farrisric/NPlib/blob/main/docs/images/logo.png?raw=true)

## Table of contents

- [NanoParticleLibrary (NPL)](#nanoparticlelibrary-npl)
  - [Table of contents](#table-of-contents)
  - [About NPL](#about-npl)
  - [Documentation](#documentation)
  - [Installation](#installation)
    - [Requirements](#requirements)
    - [Installation from PyPI](#installation-from-pypi)
    - [Installation from source](#installation-from-source)
  - [Usage](#usage)
    - [Local Optimization](#local-optimization)
    - [Basin-Hopping](#basin-hopping)
  - [Examples](#examples)
  - [Development](#development)
  - [References](#references)
  - [Contact](#contact)
  - [License](#license)

## About NPL

NPL is a Python library for the simulation and structural optimization of nanoparticles, specifically tailored for bimetallic nanoparticles. Built on the robust ASE (Atomic Simulation Environment), it enables users to easily set up and analyze complex nanoparticle structures across a range of chemical compositions and structures. NPL provides high-level abstractions, making it accessible for both beginners and experienced researchers aiming to perform detailed nanoparticle simulations.

## Documentation

A partial documentation is available at: https://nplib.readthedocs.io/en/latest/

## Installation

### Requirements

- Python 3.9+
- Atomic Simulation Environment (ASE) >= 3.21
- scikit-learn
- sortedcontainers

### Installation from PyPI

You can install NPL with pip:

```sh
pip install npl
```

or from github:

```sh
git clone https://github.com/farrisric/NPlib
pip install ./NPlib
```

## Examples

### Monte Carlo Run Example

Here is an example of how to perform a Monte Carlo run using NPL:

```python

from npl.descriptors import ExtendedTopologicalFeaturesClassifier
from npl.calculators import TOPCalculator
from npl.core import Nanoparticle
from npl.monte_carlo import run_monte_carlo
from npl.visualize import plot_parted_particle

energy_calculator = TOPCalculator('ETOP', stoichiometry='Pt151Cu50',
                     feature_classifier=ExtendedTopologicalFeaturesClassifier)

feature_classifier = calc.get_feature_classifier()

beta = 250
max_steps = 10000

start_particle = Nanoparticle()
start_particle.truncated_octahedron(7, 2, {'Pt': 151, 'Cu': 50})
best_particle, accepted_energies = run_monte_carlo(beta, max_steps,
                                                    start_particle,
                                                    energy_calculator,
                                                    feature_classifier)

plot_parted_particle(best_particle)
```

![Tutorial Image](https://github.com/farrisric/NPlib/blob/main/docs/images/tutorial4_image1.png?raw=true)

This example load pre-trained Topological coefficients, initializes a truncated octahedral Pt151Cu50 nanoparticle, sets up a Monte Carlo simulation at beta 250 for 10000 steps, runs the simulation, and then prints the optimized positions of the particle.

## References

If you use this code, please cite our papers:

```bibtex
@neuman{10.1063/5.0214377,
    author = {Felix Neumann  and Johannes T Margraf and Karsten Reuter and Albert Bruix},
    title = "{Interplay between shape and composition in bimetallic nanoparticles
    revealed by an efficient optimal-exchange optimization algorithm}",
    archivePrefix = {ChemRxiv},
    doi = {10.26434/chemrxiv-2021-26ztp},
}

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

@farris{10.1063/5.0214377,
    author = {Farris, Riccardo and Neyman, Konstantin M. and Bruix, Albert},
    title = "{Determining the chemical ordering in nanoalloys by considering atomic coordination types}",
    journal = {The Journal of Chemical Physics},
    volume = {161},
    number = {13},
    pages = {134114},
    year = {2024},
    issn = {0021-9606},
    doi = {10.1063/5.0214377}
}
```
