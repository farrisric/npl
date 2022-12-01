from numba import njit
import numpy as np

from npl.descriptors import SiteFinder
from npl.core import Nanoparticle

class SiteDescriptor(SiteFinder):
    """Class for the calculation of adsoption site descriptors as in Mie et al. 2019.
    Each site has its own set of descriptors.
    
    Parameters
    ----------
    system:
        identifier of the descriptor
    """

    def __init__(self, particle):
        super().__init__(particle)

    def create(self, particle):
        return super().create(particle)

    
    

    