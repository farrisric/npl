import pytest
from npl.core import Nanoparticle


@pytest.fixture
def ptcu_particle():
    """A 201-atom truncated-octahedral Pt151Cu50 nanoparticle."""
    particle = Nanoparticle()
    particle.truncated_octahedron(7, 2, {"Pt": 151, "Cu": 50})
    return particle
