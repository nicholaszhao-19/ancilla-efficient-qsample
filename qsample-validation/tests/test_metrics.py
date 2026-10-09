import numpy as np
import pytest

from qsample_validation.metrics import (
    clean_joint_state,
    clean_subspace_leakage,
    joint_trace_distance,
    normalized_state,
)


def test_small_joint_error_and_leakage_avoid_subtractive_cancellation() -> None:
    probabilities = np.array([1.0, 0.0])
    state = clean_joint_state(probabilities)
    leaked_amplitude = 1.0e-9
    state[1, 0, 0] = leaked_amplitude

    assert joint_trace_distance(state, probabilities) == pytest.approx(
        leaked_amplitude, rel=1.0e-12
    )
    assert clean_subspace_leakage(state) == pytest.approx(
        leaked_amplitude**2, rel=1.0e-12
    )


def test_state_normalization_is_scale_robust() -> None:
    state = np.array([1.0, 0.1j], dtype=np.complex128)
    expected = normalized_state(state)
    np.testing.assert_allclose(normalized_state(state * 1.0e-200), expected)
    np.testing.assert_allclose(normalized_state(state * 1.0e200), expected)
