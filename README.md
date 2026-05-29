<!-- [![GitHub release](https://img.shields.io/github/release/yourusername/npl.svg)](https://GitHub.com/yourusername/npl/releases/) -->

[![License](https://img.shields.io/badge/License-MIT%202.0-blue.svg)](https://opensource.org/licenses/mit)
[![GitHub issues](https://img.shields.io/github/issues/farrisric/npl.svg)](https://GitHub.com/farrisric/npl/issues)
[![Documentation Status](https://readthedocs.org/projects/nplib/badge/)](https://nplib.readthedocs.io/en/latest/index.html)
[![DOI](https://sandbox.zenodo.org/badge/549641845.svg)](https://handle.stage.datacite.org/10.5072/zenodo.262134)

# <span style="font-size:larger;">NanoParticleLibrary (NPL)</span>

![NPL Logo](https://github.com/farrisric/npl/blob/main/docs/images/logo.png?raw=true)

## Table of contents

- [NanoParticleLibrary (NPL)](#nanoparticlelibrary-npl)
  - [About NPL](#about-npl)
  - [Documentation](#documentation)
  - [Installation](#installation)
    - [Requirements](#requirements)
    - [Installation from PyPI](#installation-from-pypi)
    - [Installation from source](#installation-from-source)
  - [Examples](#examples)
  - [References](#references)
  - [Contact](#contact)
  - [License](#license)

## About NPL

NPL is a Python library for the simulation and structural optimization of nanoparticles, specifically tailored for bimetallic nanoparticles. Built on the robust ASE (Atomic Simulation Environment), it enables users to easily set up and analyze complex nanoparticle structures across a range of chemical compositions and structures. NPL provides high-level abstractions, making it accessible for both beginners and experienced researchers aiming to perform detailed nanoparticle simulations.

## Documentation

A partial documentation is available at: https://nplib.readthedocs.io/en/latest/

## Installation

### Requirements

- Python 3.10+
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
git clone https://github.com/farrisric/npl
pip install ./npl
```

## Examples

### TOP Energy Evaluation Example

This example demonstrates how to evaluate the topological (TOP) energy of a bimetallic nanoparticle using pretrained coefficients shipped with NPL. A 201-atom truncated-octahedral Pt151Cu50 particle is constructed, its extended topological feature vector is computed, and the energy is evaluated using a `TOPCalculator` loaded from pretrained linear-model coefficients. The result is a fast, DFT-quality energy estimate without requiring an external calculator.

```python
from npl.core import Nanoparticle
from npl.calculators import TOPCalculator
from npl.descriptors import ExtendedTopologicalFeaturesClassifier

calculator = TOPCalculator(
    "ETOP",
    stoichiometry="Pt151Cu50",
    feature_classifier=ExtendedTopologicalFeaturesClassifier,
)
feature_classifier = calculator.get_feature_classifier()

particle = Nanoparticle()
particle.truncated_octahedron(7, 2, {"Pt": 151, "Cu": 50})

feature_classifier.compute_feature_vector(particle)
energy = calculator.compute_energy(particle)

print(f"TOP energy of Pt151Cu50: {float(energy):.4f}")
```

> **Looking for Monte Carlo / grand-canonical sampling?** That functionality now lives in the companion package [mcpy](https://github.com/farrisric/mcpy).

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

## Contact

For any questions or issues, please contact:

- Riccardo Farris: [rfarris@ub.edu](mailto:rfarris@ub.edu)
- GitHub Issues: [npl Issues](https://github.com/farrisric/npl/issues)

## License

This project is licensed under the MIT License - see the [LICENSE](https://opensource.org/licenses/MIT) file for details.
