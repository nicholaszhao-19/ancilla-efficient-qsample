import numpy as np

from qsample_validation.gqsp_filter import (
    average_polynomial,
    design_filter,
    evaluate_filter,
    scalar_completion,
)


def test_average_power_filter_and_completion() -> None:
    gap = 0.2
    eta = 1.0e-5
    spec = design_filter(gap, eta)
    assert abs(evaluate_filter(1.0 + 0.0j, spec) - 1.0) < 1.0e-12

    phases = np.linspace(gap, np.pi, 20_000)
    phases = np.concatenate((-phases, phases))
    values = evaluate_filter(np.exp(1.0j * phases), spec)
    assert np.max(np.abs(values)) <= eta * (1.0 + 1.0e-10)

    for phase in np.linspace(-np.pi, np.pi, 31):
        block = scalar_completion(np.exp(1.0j * phase), spec)
        np.testing.assert_allclose(
            block.conj().T @ block, np.eye(2), atol=1.0e-12
        )
        assert abs(block[0, 0] - evaluate_filter(np.exp(1.0j * phase), spec)) < 1e-12


def test_filter_is_stable_at_tiny_representable_gaps() -> None:
    gap = 5.0e-14
    eta = 1.0e-6
    spec = design_filter(gap, eta)
    value_at_gap = abs(evaluate_filter(np.exp(1.0j * gap), spec))
    assert value_at_gap <= spec.certified_eta * (1.0 + 1.0e-10)
    assert average_polynomial(0.0 + 0.0j, 7) == 1.0 / 7.0

    smallest_eta = float(np.nextafter(0.0, 1.0))
    extreme = design_filter(0.2, smallest_eta)
    assert extreme.n == int(np.ceil(-np.log(smallest_eta)))
