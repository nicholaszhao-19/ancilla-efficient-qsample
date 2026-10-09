import numpy as np

from qsample_validation.ising import IsingModel


def test_open_ladder_geometry_and_energy() -> None:
    model = IsingModel(2, 4)
    assert model.n_states == 256
    assert len(model.edges) == 10
    energies = model.energies()
    assert energies.min() == 0.0
    assert energies.max() == 10.0
    assert np.count_nonzero(energies == 0.0) == 2
    unique, counts = np.unique(energies, return_counts=True)
    np.testing.assert_array_equal(unique, [0, 2, 3, 4, 5, 6, 7, 8, 10])
    np.testing.assert_array_equal(counts, [2, 14, 24, 48, 80, 48, 24, 14, 2])


def test_configuration_index_round_trip() -> None:
    model = IsingModel(2, 3)
    configurations = model.configurations()
    np.testing.assert_array_equal(
        model.index_from_spins(configurations), np.arange(model.n_states)
    )
