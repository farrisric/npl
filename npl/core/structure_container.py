
import numpy as np
from ase import Atoms

from npl.descriptors import BaseDescriptor

class StructureContainer:
    """
    Class for structure storing to train surrogate energy models"""    

    def __init__(self, descriptor : BaseDescriptor):
        self._structure_list = []
        self._descriptor = descriptor

    def __len__(self) -> int:
        return len(self._structure_list)

    def __getitem__(self, ind: int):
        return self._structure_list[ind]

    def add_structure(self, structure: Atoms, user_tag: str = None,
                      properties: dict = None, allow_duplicate: bool = True,
                      sanity_check: bool = True):
        """
        Adds a structure to the structure container.

        Parameters
        ----------
        structure
            the atomic structure to be added
        user_tag
            custom user tag to label structure
        properties
            scalar properties. If properties are not specified the structure
            object will be checked for an attached ASE calculator object
            with a calculated potential energy
        allow_duplicate
             whether or not to add the structure if there already exists a
             structure with identical cluster-vector
        sanity_check
            whether or not to carry out a sanity check before adding the
            structure. This includes checking occupations and volume.
        """

        # structure must have a proper format and label
        if not isinstance(structure, Atoms):
            raise TypeError('structure must be an ASE Atoms object not {}'.format(type(structure)))

        if user_tag is not None:
            if not isinstance(user_tag, str):
                raise TypeError('user_tag must be a string not {}.'.format(type(user_tag)))

        if sanity_check:
            pass #write something for sanity check

        # check for properties in attached calculator
        if properties is None:
            properties = {}
            if structure.calc is not None:
                if not structure.calc.calculation_required(structure, ['energy']):
                    energy = structure.get_potential_energy()
                    properties['energy'] = energy

        # check if there exist structures with identical cluster vectors
        structure_copy = structure.copy()
        feature_vector = self._descriptor.get_feature_vector(structure_copy)
        if not allow_duplicate:
            for i, fs in enumerate(self):
                if np.allclose(feature_vector, fs.feature_vector):
                    msg = '{} and {} have identical feature vectors'.format(
                        user_tag if user_tag is not None else 'Input structure',
                        fs.user_tag if fs.user_tag != 'None' else 'structure')
                    msg += ' at index {}'.format(i)
                    raise ValueError(msg)

        # add structure
        structure = FitStructure(structure_copy, user_tag, feature_vector, properties)
        self._structure_list.append(structure)


class FitStructure:
    """
    This class holds a supercell along with its properties and cluster
    vector.

    Attributes
    ----------
    structure : Atoms
        supercell structure
    user_tag : str
        custom user tag
    cvs : np.ndarray
        calculated cluster vector for actual structure
    properties : dict
        dictionary of properties
    """

    def __init__(self, structure: Atoms, user_tag: str,
                 feature_vector: np.ndarray, properties: dict = {}):
        self._structure = structure
        self._user_tag = user_tag
        self._feature_vector = feature_vector
        self.properties = properties

    @property
    def cluster_vector(self) -> np.ndarray:
        """calculated feature vector"""
        return self._feature_vector

    @property
    def structure(self) -> Atoms:
        """atomic structure"""
        return self._structure

    @property
    def user_tag(self) -> str:
        """structure label"""
        return str(self._user_tag)

    def __getattr__(self, key):
        """ Accesses properties if possible and returns value. """
        if key not in self.properties.keys():
            return super().__getattribute__(key)
        return self.properties[key]

    def __len__(self) -> int:
        """ Number of sites in the structure. """
        return len(self._structure)
