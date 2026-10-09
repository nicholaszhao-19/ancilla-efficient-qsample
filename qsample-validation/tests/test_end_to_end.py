import pytest

import qsample_validation.end_to_end as end_to_end
from qsample_validation.end_to_end import run_end_to_end
from qsample_validation.ising import IsingModel
from qsample_validation.markov import build_chain


def test_complete_two_by_four_schedule_meets_requested_error() -> None:
    model = IsingModel(2, 4)
    chains = [build_chain(model, beta) for beta in (0.0, 0.3, 0.6, 0.9, 1.2)]
    epsilon = 0.05
    result = run_end_to_end(chains, epsilon)
    assert result.joint_error <= epsilon + 1.0e-8
    assert result.reduced_error <= epsilon + 1.0e-8
    assert result.total_variation <= epsilon + 1.0e-8
    assert result.norm_residual < 1.0e-10
    assert result.reduced_error <= result.joint_error + 1.0e-10
    assert result.total_variation <= result.reduced_error + 1.0e-10
    for stage in result.stages:
        expected_eta = epsilon / (4.0 * 4 * stage.phase_weight)
        assert stage.eta == pytest.approx(expected_eta)


def test_uniform_filter_allocation_matches_corrected_proceedings() -> None:
    model = IsingModel(2, 4)
    chains = [build_chain(model, beta) for beta in (0.0, 0.3, 0.6, 0.9, 1.2)]
    epsilon = 1.0e-3
    result = run_end_to_end(chains, epsilon, filter_allocation="uniform")
    assert result.joint_error <= epsilon + 1.0e-10
    assert result.reduced_error <= result.joint_error + 1.0e-10
    assert result.total_variation <= result.reduced_error + 1.0e-10
    for stage in result.stages:
        expected = epsilon / (4.0 * 4 * (stage.sequence.length - 1))
        assert stage.filter_allocation == "uniform"
        assert stage.eta == pytest.approx(expected)


def test_paper_benchmark_headline_metrics_and_resources(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    model = IsingModel(2, 4)
    chains = [build_chain(model, beta) for beta in (0.0, 0.3, 0.6, 0.9, 1.2)]
    compiled_calls: list[tuple[float, float]] = []
    compiled_phase = end_to_end.apply_compiled_selective_phase

    def traced_compiled_phase(state, walk, spec, phase):
        compiled_calls.append((walk.chain.beta, float(phase)))
        return compiled_phase(state, walk, spec, phase)

    monkeypatch.setattr(
        end_to_end,
        "apply_compiled_selective_phase",
        traced_compiled_phase,
    )
    result = run_end_to_end(chains, 1.0e-3)

    assert result.joint_error == pytest.approx(1.08167838e-4, rel=5.0e-8)
    assert result.ideal_joint_error == pytest.approx(1.08167837462e-4, rel=5.0e-10)
    assert result.compiled_to_ideal_joint_distance == pytest.approx(1.08475099934e-8, rel=5.0e-8)
    assert result.compiled_to_ideal_joint_distance > 0.0
    assert result.joint_error <= (
        result.ideal_joint_error + result.compiled_to_ideal_joint_distance + 1.0e-12
    )
    assert result.reduced_error <= result.joint_error + 1.0e-12
    assert result.total_variation <= result.reduced_error + 1.0e-12
    assert 0.0 < result.leakage < 1.0e-14
    assert result.norm_residual < 1.0e-10
    assert sum(stage.resources.controlled_walk_calls for stage in result.stages) == 10516
    assert sum(stage.resources.wrapper_oracle_calls for stage in result.stages) == 32
    assert sum(stage.resources.transition_oracle_calls for stage in result.stages) == 21064
    assert [stage.phase_sensitive_compiler_bound for stage in result.stages] == pytest.approx(
        [1.66073e-7, 2.47241e-7, 2.33881e-7, 1.14022e-7],
        rel=5.0e-5,
    )
    assert result.compiled_to_ideal_joint_distance <= sum(
        stage.phase_sensitive_compiler_bound for stage in result.stages
    )

    expected_calls: list[tuple[float, float]] = []
    for stage in result.stages:
        word = []
        for alpha, beta in zip(stage.sequence.alpha, stage.sequence.beta, strict=True):
            word.extend(
                [
                    (stage.beta_target, float(beta)),
                    (stage.beta_source, float(-alpha)),
                ]
            )
        expected_calls.extend(word)
        expected_calls.extend(word)
    assert len(compiled_calls) == 2 * sum(stage.sequence.n_phases for stage in result.stages)
    assert [item[0] for item in compiled_calls] == pytest.approx(
        [item[0] for item in expected_calls]
    )
    assert [item[1] for item in compiled_calls] == pytest.approx(
        [item[1] for item in expected_calls]
    )
