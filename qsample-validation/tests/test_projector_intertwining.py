from qsample_validation.ising import IsingModel
from qsample_validation.markov import build_chain
from qsample_validation.projectors import (
    projector_intertwining_operator_norm,
    projector_intertwining_probe_residual,
)
from qsample_validation.szegedy import SpectralWalk


def test_projector_intertwining_on_three_ladders() -> None:
    for cols in (2, 3, 4):
        walk = SpectralWalk(build_chain(IsingModel(2, cols), 0.7))
        assert projector_intertwining_probe_residual(walk, seed=cols) < 1.0e-10
        assert projector_intertwining_operator_norm(walk) < 1.0e-10
