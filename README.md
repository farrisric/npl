# <span style="font-size:larger;">NanoParticleLibrary (NPL)</span>

<!-- [![GitHub release](https://img.shields.io/github/release/yourusername/npl.svg)](https://GitHub.com/yourusername/npl/releases/) -->
[![License](https://img.shields.io/badge/License-MIT%202.0-blue.svg)](https://opensource.org/licenses/mit)
[![GitHub issues](https://img.shields.io/github/issues/farrisric/nplib.svg)](https://GitHub.com/farrisric/NPlib/issues)
[![Documentation Status](https://readthedocs.org/projects/nplib/badge/)](https://nplib.readthedocs.io/en/latest/modules.html)

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

NPL (NanoParticleLibrary) is a comprehensive wrapper around the popular ASE (Atomic Simulation Environment) library, designed to facilitate the manipulation and optimization of nanoparticles. It is particularly tailored for optimizing the chemical ordering in bimetallic nanoparticles.

## Documentation

A partial documentation is available at: https://npl-docs.readthedocs.io

## Installation

### Requirements

- Python 3.7+
- Atomic Simulation Environment (ASE) >= 3.21
- scikit-learn
- sortedcontainers

### Installation from PyPI

This is the recommended way to install NPL.

```sh
git clone https://github.com/farrisric/NPlib
pip install ./NPlib
```
## References

If you use this code, please cite our papers:

```bibtex
@farris2024{10.1063/5.0214377,
    author = {Farris, Riccardo and Neyman, Konstantin M. and Bruix, Albert},
    title = "{Determining the chemical ordering in nanoalloys by considering atomic coordination types}",
    journal = {The Journal of Chemical Physics},
    volume = {161},
    number = {13},
    pages = {134114},
    year = {2024},
    issn = {0021-9606},
    doi = {10.1063/5.0214377},
}
```