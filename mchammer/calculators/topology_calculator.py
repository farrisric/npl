
from ase import Atoms
from npl.descriptors.topologies import Topologies
from mchammer.calculators.base_calculator import BaseCalculator

class TopologicalCalculator(BaseCalculator):
    """A Topological Calculator object that enables the 
    calculation of the nanoparticle energy based on the
    chemical ordering
    
    Parameters
    ----------
    structure : ase.Atoms
        structure for which to set up the calculator
    topology : 
    """
    def __init__(self,
                 structure: Atoms, topology: Topologies,
                 name: str = 'Topological Calculator') -> None:
        super().__init__(name=name)

        structure_cpy = structure.copy()
        self._topology = topology


    @property
    def topology(self) -> Topologies:
        """ topology  """
        return self._cluster_expansion