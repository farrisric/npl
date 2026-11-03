# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
import tomllib
sys.path.insert(0, os.path.abspath('..'))

project = 'Nanoparticle Library'
copyright = '2024, Riccardo Farris'
author = 'Riccardo Farris'
with open(os.path.join(os.path.dirname(__file__), '..', 'pyproject.toml'), 'rb') as fh:
    release = tomllib.load(fh)['project']['version']

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc', 'sphinx.ext.napoleon', 'sphinx.ext.todo', 'sphinx.ext.viewcode']

# package __init__ files re-export their submodules' classes; document each class once
autodoc_default_options = {'ignore-module-all': True}
# class docstrings list their methods; napoleon would otherwise emit them a second time
napoleon_custom_sections = [('Methods', 'params_style')]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = []
