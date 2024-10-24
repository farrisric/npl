# <span style="font-size:larger;">NanoParticleLibrary (NPL)</span>

[![GitHub release](https://img.shields.io/github/release/yourusername/npl.svg)](https://GitHub.com/yourusername/npl/releases/)
[![License](https://img.shields.io/badge/License-MIT%202.0-blue.svg)](https://opensource.org/licenses/mit)
[![GitHub issues](https://img.shields.io/github/issues/yourusername/npl.svg)](https://GitHub.com/yourusername/npl/issues/)
[![Documentation Status](https://readthedocs.org/projects/npl/badge/)](https://nplib.readthedocs.io/en/latest/modules.html)

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
pip install --upgrade pip
pip install nanoparticlelibrary