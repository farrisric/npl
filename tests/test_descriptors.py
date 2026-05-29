import numpy as np
from npl.core import Nanoparticle
from npl.descriptors import (
    ExtendedTopologicalFeaturesClassifier,
    NeighborCountingEnvironmentCalculator,
)


def test_neighbor_counting_environment(ptcu_particle):
    calc = NeighborCountingEnvironmentCalculator(["Pt", "Cu"])
    env = calc.predict_local_environment(ptcu_particle, 0)
    assert len(env) == 2
    assert env.sum() == len(ptcu_particle.neighbor_list[0])


def test_extended_topological_feature_vector(ptcu_particle):
    classifier = ExtendedTopologicalFeaturesClassifier(["Pt", "Cu"])
    classifier.compute_feature_vector(ptcu_particle)
    key = classifier.get_feature_key()
    fv = ptcu_particle.get_feature_vector(key)
    assert fv is not None
    assert np.asarray(fv).ndim == 1
    assert len(fv) == len(classifier.get_feature_labels())


def test_compute_atom_feature_for_symbol_matches_current_symbol():
    fc = ExtendedTopologicalFeaturesClassifier(["Au", "Cu", "Pt"])
    p = Nanoparticle()
    p.truncated_octahedron(7, 2, {"Pt": 101, "Au": 50, "Cu": 50})
    fc.compute_feature_vector(p)  # populates sublayer_indices
    idx = p.get_indices_by_symbol("Pt")[0]
    same = fc.compute_atom_feature_for_symbol(p, idx, "Pt")
    assert np.array_equal(same, fc.compute_atom_feature(p, idx))
    # asking for a different symbol changes the feature
    other = fc.compute_atom_feature_for_symbol(p, idx, "Au")
    assert not np.array_equal(same, other)
