def test_truncated_octahedron_atom_count(ptcu_particle):
    assert ptcu_particle.atoms.get_n_atoms() == 201


def test_truncated_octahedron_stoichiometry(ptcu_particle):
    stoich = ptcu_particle.get_stoichiometry()
    assert stoich["Pt"] == 151
    assert stoich["Cu"] == 50


def test_neighbor_list_is_built(ptcu_particle):
    indices = list(ptcu_particle.get_indices())
    assert len(indices) == 201
    assert all(len(ptcu_particle.neighbor_list[i]) > 0 for i in indices)


def test_geometrical_data_has_positions_and_neighbors(ptcu_particle):
    data = ptcu_particle.get_geometrical_data()
    assert "positions" in data and "neighbor_list" in data
    assert len(data["positions"]) == 201
