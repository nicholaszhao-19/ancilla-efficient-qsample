import numpy as np

from qsample_validation.ising import IsingModel
from qsample_validation.markov import build_chain, squared_qsample_overlap


def test_exact_benchmark_values_and_residuals() -> None:
    model = IsingModel(2, 4)
    beta = [0.0, 0.3, 0.6, 0.9, 1.2]
    chains = [build_chain(model, value) for value in beta]
    expected_gaps = np.array(
        [0.062500, 0.039040, 0.021008, 0.009999, 0.004360]
    )
    expected_phase_gaps = np.array(
        [0.355421, 0.280346, 0.205340, 0.141532, 0.093413]
    )
    expected_overlaps = np.array([0.944415, 0.939133, 0.933454, 0.936311])
    np.testing.assert_allclose(
        [chain.spectral_gap for chain in chains], expected_gaps, atol=5.0e-7
    )
    np.testing.assert_allclose(
        [chain.walk_phase_gap for chain in chains],
        expected_phase_gaps,
        atol=5.0e-7,
    )
    np.testing.assert_allclose(
        [
            squared_qsample_overlap(
                chains[index].probabilities, chains[index + 1].probabilities
            )
            for index in range(4)
        ],
        expected_overlaps,
        atol=5.0e-7,
    )
    for chain in chains:
        assert chain.row_residual < 1.0e-12
        assert chain.stationarity_residual < 1.0e-12
        assert chain.detailed_balance_residual < 1.0e-12

