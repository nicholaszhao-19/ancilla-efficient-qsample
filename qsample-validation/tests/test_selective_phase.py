import numpy as np

from qsample_validation.gqsp_filter import design_filter
from qsample_validation.ising import IsingModel
from qsample_validation.markov import build_chain
from qsample_validation.metrics import clean_joint_state
from qsample_validation.selective_phase import (
    apply_compiled_selective_phase,
    apply_filter_completion,
    exact_restricted_error,
)
from qsample_validation.szegedy import SpectralWalk


def test_complete_compiler_restricted_operator_error() -> None:
    chain = build_chain(IsingModel(2, 2), 0.8)
    walk = SpectralWalk(chain)
    eta = 1.0e-4
    spec = design_filter(chain.walk_phase_gap, eta)
    phase = np.pi
    stationary = np.sqrt(chain.probabilities).astype(np.complex128)
    columns = []
    for index in range(walk.n_states):
        logical = np.zeros(walk.n_states, dtype=np.complex128)
        logical[index] = 1.0
        clean = np.zeros((2, walk.n_states, walk.n_states), dtype=np.complex128)
        clean[0, :, 0] = logical
        actual = apply_compiled_selective_phase(clean, walk, spec, phase)
        ideal_logical = logical + (
            np.exp(1.0j * phase) - 1.0
        ) * stationary * np.vdot(stationary, logical)
        ideal = clean_joint_state(np.abs(ideal_logical) ** 2)
        ideal[0, :, 0] = ideal_logical
        columns.append((actual - ideal).reshape(-1))
    operator = np.stack(columns, axis=1)
    direct_error = float(np.linalg.svd(operator, compute_uv=False)[0])
    formula_error = exact_restricted_error(walk, spec, phase)
    np.testing.assert_allclose(direct_error, formula_error, atol=1.0e-10)
    assert direct_error <= 2.0 * eta + 1.0e-10


def test_spectral_completion_and_compiler_are_unitary_off_clean_column() -> None:
    chain = build_chain(IsingModel(2, 2), 0.8)
    walk = SpectralWalk(chain)
    spec = design_filter(chain.walk_phase_gap, 1.0e-4)
    generator = np.random.default_rng(22)
    state = generator.normal(size=(2, walk.n_states, walk.n_states))
    state = state + 1.0j * generator.normal(size=state.shape)
    state = state.astype(np.complex128)
    state /= np.linalg.norm(state)

    computed = apply_filter_completion(state, walk, spec)
    recovered = apply_filter_completion(computed, walk, spec, adjoint=True)
    np.testing.assert_allclose(recovered, state, atol=3.0e-12)

    phased = apply_compiled_selective_phase(state, walk, spec, 0.73)
    unphased = apply_compiled_selective_phase(phased, walk, spec, -0.73)
    np.testing.assert_allclose(unphased, state, atol=5.0e-12)
