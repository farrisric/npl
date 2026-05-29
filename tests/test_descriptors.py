import numpy as np
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
