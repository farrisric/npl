---
title: 'NanoParticleLibrary: topological descriptors and chemical-ordering optimization for bimetallic nanoparticles'
tags:
  - Python
  - computational chemistry
  - materials science
  - nanoparticles
  - nanoalloys
  - chemical ordering
  - heterogeneous catalysis
authors:
  - name: Riccardo Farris
    orcid: 0000-0002-8157-6786
    affiliation: "1, 2"
  - name: Albert Bruix
    orcid: 0000-0003-2585-5542
    affiliation: 1
affiliations:
  - name: Institut de Química Teòrica i Computacional (IQTCUB), Universitat de Barcelona, Barcelona, Spain
    index: 1
  - name: LEITAT Technological Center, Terrassa, Spain
    index: 2
date: 29 May 2026
bibliography: paper.bib
---

# Summary

`NanoParticleLibrary` (`npl`) is a Python package for building bimetallic
nanoparticles, describing them through coordination-based topological
descriptors, and optimizing their chemical ordering. It is built on the Atomic
Simulation Environment (ASE) [@ase2017], representing each particle through a
neighbor list and per-atom local environments from which feature vectors are
computed. These feature vectors feed a linear topological energy model whose
coefficients are fitted with scikit-learn [@scikit-learn] regressors (Bayesian
ridge regression) against reference energies, or loaded from pretrained
parameters shipped with the package. Given an energy model, `npl` searches for
low-energy chemical orderings with guided local optimization, basin hopping,
and genetic-algorithm strategies. The topological descriptors and the
optimal-exchange ideas behind these methods were introduced in
[@neumann2021; @farris2024ordering] and applied to catalytic nanoparticles in
[@farris2024zr].

# Statement of need

The properties of bimetallic nanoparticles depend not only on their size and
shape but strongly on their chemical ordering — how the two elements are
arranged among the available sites. Exhaustively enumerating orderings is
combinatorially impossible, and evaluating each candidate with density
functional theory is prohibitive. `npl` addresses this by combining cheap,
physically motivated topological descriptors with fitted energy models and
dedicated ordering optimizers, so that researchers can predict and explore the
most stable arrangements of nanoalloys at low cost. It is aimed at
computational chemists and materials scientists studying nanoalloys for
heterogeneous catalysis, where chemical ordering controls activity and
stability.

`npl` focuses on structure and chemical-ordering optimization with a fixed
number of atoms. Sampling at controlled chemical potential — grand canonical
Monte Carlo with machine-learning potentials — is provided by the companion
package `mcpy` [@mcpy]; the two packages share an ASE-native design and are
intended to be used together.

# Software design

`npl` is organized around four extensible abstractions connected by string
keys on a central particle object. A `Nanoparticle` holds atoms, a neighbor
list, and dictionaries — keyed by name — of energies, feature vectors, and
local environments, so that several descriptors and energy models can coexist
on one particle and geometry can be reused across different orderings.
*Descriptors* compute global feature vectors (topological and
coordination-based classifiers) or per-atom local environments. *Calculators*
implement a common `compute_energy(particle)` interface and write into the
particle's energy dictionary; they include an ASE EMT calculator, scikit-learn
regressors over feature vectors, and the linear topological model with
pretrained coefficients. *Optimizers* search chemical-ordering space using
guided exchanges, basin hopping, and genetic algorithms. This separation lets
users swap energy models and search strategies without touching the particle
representation.

# Acknowledgements

The authors received no specific funding for this work.

# AI usage disclosure

Claude (Opus) assisted with code modernization, dead-code removal, test
scaffolding, packaging configuration, and editing of this manuscript and the
documentation. All AI-assisted outputs were reviewed, edited, and validated by
the human authors, who designed the code architecture and take full
responsibility for the accuracy and originality of the software and all
submitted materials.

# References
