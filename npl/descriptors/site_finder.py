from numba import njit
import numpy as np

from npl.descriptors import BaseDescriptor
from npl.core import Nanoparticle

class SiteFinder(BaseDescriptor):
    """Class for the calculation of adsoption site descriptors as in Mie et al. 2019.
    Each site has its own set of descriptors.
    
    Parameters
    ----------
    system:
        identifier of the descriptor
    """
    
    def __init__(self, particle):
        self.particle = particle
        self.surface_neighbors_dict = None
        self.singlets = []
        self.pairs = []
        self.triplets = []
        self.quadruplets = []
        
        super().__init__(name='ADS')
        self.get_site_singlets()
        self.constructur_surface_neighbors()
        self.get_site_pairs()
        self.get_site_triplets()
        self.get_site_quadruplets()

        del self.particle

    def create(self, particle):
        pass 

    def constructur_surface_neighbors(self):
        self.surface_neighbors_dict = {x : [] for x in self.singlets}
        
        for surface_idx in self.singlets:
            for neigh in self.particle.neighbor_dict[surface_idx]:
                if neigh in self.singlets:
                    self.surface_neighbors_dict[surface_idx].append(neigh)

    def get_site_singlets(self):
        coordination_list = self.particle.get_coordination_list()
        self.singlets = np.where(coordination_list<12)[0].tolist()

    def get_site_pairs(self):

        for central_i, neighbors in self.surface_neighbors_dict.items():
            for neigh_i in sorted(neighbors):
                if neigh_i > central_i:
                    self.pairs.append([central_i, neigh_i])

    def get_site_triplets(self):

        for idx1, idx2 in self.pairs:
            common_idices = np.intersect1d(self.surface_neighbors_dict[idx1], self.surface_neighbors_dict[idx2])
            for common_idx in common_idices:
                site_indices = sorted([idx1, idx2,common_idx])
                if site_indices not in self.triplets:
                    self.triplets.append(site_indices)

    def get_site_quadruplets(self):

        for i, bridge_site in enumerate(self.pairs):
            for other_bridge_site in self.pairs[i+1:]:
                if len(np.intersect1d(bridge_site,other_bridge_site)) == 0:
                    a, b = bridge_site
                    c, d = other_bridge_site 

                    check_a =  np.intersect1d(self.surface_neighbors_dict[a], other_bridge_site)
                    check_b = np.intersect1d(self.surface_neighbors_dict[b], other_bridge_site)

                    if len(check_a) == 1 and len(check_b) == 1:
                        if check_a != check_b:
                            site_indices = sorted([a,b,c,d])
                            if site_indices not in self.quadruplets:
                                self.quadruplets.append(site_indices)





    

    