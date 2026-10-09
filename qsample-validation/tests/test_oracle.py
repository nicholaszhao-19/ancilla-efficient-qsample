import numpy as np

from qsample_validation.ising import IsingModel
from qsample_validation.markov import build_chain
from qsample_validation.oracle import HouseholderOracle


def test_householder_clean_column_and_unitarity() -> None:
    chain = build_chain(IsingModel(2, 2), 0.7)
    oracle = HouseholderOracle.from_transition(chain.transition)
    clean = np.zeros((oracle.n_states, oracle.n_states), dtype=np.complex128)
    clean[:, 0] = 1.0
    prepared = oracle.apply(clean)
    np.testing.assert_allclose(prepared, np.sqrt(chain.transition), atol=1.0e-12)

    generator = np.random.default_rng(4)
    state = generator.normal(size=clean.shape) + 1.0j * generator.normal(
        size=clean.shape
    )
    recovered = oracle.apply(oracle.apply(state))
    np.testing.assert_allclose(recovered, state, atol=1.0e-12)

