"""
NanoParticleLibrary (npl) — simulate and optimize the structure and chemical
ordering of (primarily bimetallic) nanoparticles, built on top of ASE.
"""
import logging
from importlib.metadata import PackageNotFoundError, version

from .logging_config import setup_logging  # noqa: F401  opt-in, for applications

# A library must not configure logging: the application decides where records go.
logging.getLogger(__name__).addHandler(logging.NullHandler())

try:
    __version__ = version("npl")  # single source: pyproject.toml
except PackageNotFoundError:  # source checkout on sys.path, e.g. the Read the Docs build
    __version__ = "unknown"
__maintainer__ = "Riccardo Farris"
__credits__ = 'Reactive Nanoparticles Laboratory'
__copyright__ = '2024'
__email__ = 'rfarris@ub.edu'
__status__ = 'Stable'
