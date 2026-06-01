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


def test_all_symbol_features_match_per_symbol():
    # The single-neighbor-walk variant used in the ETOP hot path must be
    # bit-identical to calling compute_atom_feature_for_symbol per species.
    symbols = ["Au", "Cu", "Pt"]
    fc = ExtendedTopologicalFeaturesClassifier(symbols)
    p = Nanoparticle()
    p.truncated_octahedron(7, 2, {"Pt": 101, "Au": 50, "Cu": 50})
    fc.compute_feature_vector(p)  # populates sublayer_indices
    for idx in [0, 50, 100, 150, 200]:
        bundle = fc.compute_atom_features_for_all_symbols(p, idx)
        assert set(bundle) == set(symbols)
        for symbol in symbols:
            ref = fc.compute_atom_feature_for_symbol(p, idx, symbol)
            assert np.array_equal(bundle[symbol], ref)
