from copy import deepcopy
from pathlib import Path

import pytest

from qsample_validation.benchmark import load_config, validate_config

CONFIG = Path(__file__).resolve().parents[1] / "configs" / "ising_2x4.yaml"


def test_frozen_benchmark_configuration_is_fully_validated() -> None:
    config = load_config(CONFIG)
    validate_config(config)

    invalid_boundary = deepcopy(config)
    invalid_boundary["model"]["boundary"] = "periodic"
    with pytest.raises(ValueError, match="model.boundary"):
        validate_config(invalid_boundary)

    invalid_laziness = deepcopy(config)
    invalid_laziness["dynamics"]["laziness"] = 0.25
    with pytest.raises(ValueError, match="dynamics.laziness"):
        validate_config(invalid_laziness)

    missing_reference_epsilon = deepcopy(config)
    missing_reference_epsilon["numerics"]["sweep_epsilon"] = [0.1, 0.01]
    with pytest.raises(ValueError, match="include benchmark_epsilon"):
        validate_config(missing_reference_epsilon)
