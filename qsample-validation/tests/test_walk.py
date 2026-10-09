import numpy as np

from qsample_validation.ising import IsingModel
from qsample_validation.markov import build_chain
from qsample_validation.szegedy import SpectralWalk


def test_walk_functional_calculus_matches_matrix_free_walk() -> None:
    walk = SpectralWalk(build_chain(IsingModel(2, 2), 0.6))
    generator = np.random.default_rng(2)
    state = generator.normal(size=(walk.n_states, walk.n_states))
    state = state + 1.0j * generator.normal(size=state.shape)
    state /= np.linalg.norm(state)

    direct = walk.apply_walk(state)
    spectral = walk.apply_function(state, lambda z: z)
    np.testing.assert_allclose(spectral, direct, atol=1.0e-12)
    assert abs(np.linalg.norm(direct) - 1.0) < 1.0e-12

    recovered = walk.apply_function(direct, lambda z: 1.0 / z)
    np.testing.assert_allclose(recovered, state, atol=1.0e-12)


def test_walk_matches_dense_definition_on_small_instance() -> None:
    walk = SpectralWalk(build_chain(IsingModel(2, 2), 0.6))
    generator = np.random.default_rng(12)
    state = generator.normal(size=(walk.n_states, walk.n_states))
    state = state + 1.0j * generator.normal(size=state.shape)
    state = state.astype(np.complex128)
    state /= np.linalg.norm(state)

    isometry = np.zeros(
        (walk.n_states * walk.n_states, walk.n_states), dtype=np.float64
    )
    for x in range(walk.n_states):
        isometry[
            x * walk.n_states : (x + 1) * walk.n_states,
            x,
        ] = np.sqrt(walk.chain.transition[x])
    reflection = 2.0 * np.einsum(
        "ik,jk->ij", isometry, isometry, optimize=False
    )
    reflection -= np.eye(walk.n_states * walk.n_states)
    dense_result = np.einsum(
        "ij,j->i", reflection, state.T.reshape(-1), optimize=False
    )

    np.testing.assert_allclose(
        walk.apply_walk(state).reshape(-1), dense_result, atol=1.0e-12
    )

    repeated = state
    for _ in range(3):
        repeated = walk.apply_walk(repeated)
    spectral_cube = walk.apply_function(state, lambda z: z**3)
    np.testing.assert_allclose(spectral_cube, repeated, atol=2.0e-12)
