"""Angle-synthesis boundary for the validation package.

The current paper outputs intentionally use the deterministic
double-precision spectral-reference backend. No rotation-angle list is emitted,
and callers must not relabel the results as gate-level GQSP simulation.
"""

from __future__ import annotations


def backend_status() -> str:
    return "spectral-reference; GQSP angles not synthesized"


def synthesize_angles() -> None:
    raise NotImplementedError(
        "Angle synthesis is outside the validated backend; use spectral-reference mode."
    )
