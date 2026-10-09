import numpy as np

from qsample_validation.fpaa import (
    apply_ideal_fpaa,
    pure_trace_distance,
    yoder_sequence,
)


def test_yoder_phase_order_and_signs() -> None:
    target = np.array([1.0, 0.0], dtype=np.complex128)
    for overlap in (1.0e-3, 1.0e-2, 0.1, 0.3, 0.8):
        source = np.array(
            [np.sqrt(overlap), np.sqrt(1.0 - overlap)], dtype=np.complex128
        )
        for target_error in (0.1, 0.03):
            sequence = yoder_sequence(overlap, target_error)
            output = apply_ideal_fpaa(source, target, sequence)
            assert pure_trace_distance(output, target) <= target_error + 1.0e-10
            assert 0.0 < sequence.phase_weight <= sequence.length - 1
