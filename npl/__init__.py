"""
NanoParticleLibrary (npl) — simulate and optimize the structure and chemical
ordering of (primarily bimetallic) nanoparticles, built on top of ASE.
"""
from importlib.metadata import version

from .logging_config import setup_logging

setup_logging()

__version__ = version("npl")  # single source: pyproject.toml
__maintainer__ = "Riccardo Farris"
__credits__ = 'Reactive Nanoparticles Laboratory'
__copyright__ = '2024'
__email__ = 'rfarris@ub.edu'
__status__ = 'Stable'
