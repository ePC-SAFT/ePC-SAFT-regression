from __future__ import annotations

import math

from epcsaft import native_sdk
import epcsaft_regression
import pytest
from epcsaft_regression import _native
from epcsaft_regression.workflow import (
    _fixed_water_factor_model,
    _water_factor_native_payload,
)


def test_water_factor_contract_is_one_parameter_over_all_nabr_rows() -> None:
    specification = epcsaft_regression.FIGIEL_WATER_SOLVATION_FACTOR_V1

    assert specification.specification_id == (
        "figiel-2025-water-solvation-factor-nabr-v1"
    )
    assert specification.source_validation_commit == (
        "8944d34f7002cda1bb8760e606cc1f11696f58cd"
    )
    assert specification.fixed_born_evidence_subject_sha256 == (
        "55ea2cd69af62c45b26179cfab6939760de23058b5a7e8c880a79f67faa417ed"
    )
    assert specification.fixed_born_diameters_angstrom == (
        2.7888130173797934,
        3.4524616464076425,
        4.147266741279482,
        4.101505615791675,
        4.476998527506598,
    )
    assert specification.fixed_aqueous_kij == (
        -0.4,
        -0.3,
        -0.1,
        -0.3,
        -0.3,
        0.8,
        0.8,
        0.0,
        0.5,
        0.65,
        -0.35,
    )
    assert specification.expected_provider_fingerprint == (
        "sha256:80ec4bde74ef7af307e44098e2e495e20a5be57d7cb23452ff321cb644816783"
    )
    assert len(specification.observations) == 21
    assert {row.salt for row in specification.observations} == {"NaBr"}
    assert tuple(row.molality_mol_per_kg for row in specification.observations) == (
        0.001,
        0.002,
        0.005,
        0.01,
        0.02,
        0.05,
        0.1,
        0.6,
        0.7,
        0.8,
        0.9,
        1.0,
        1.2,
        1.4,
        3.0,
        3.5,
        4.0,
        4.5,
        5.0,
        5.5,
        6.0,
    )
    assert specification.parameter_name == "water_solvation_factor"
    assert specification.parameter_bounds == (1.0, 2.0)
    assert specification.starts == (1.2, 1.8)
    assert specification.max_num_iterations == 500
    assert specification.start_wall_time_max_seconds == 180.0
    assert specification.function_tolerance == 1.0e-10
    assert specification.gradient_tolerance == 1.0e-10
    assert specification.parameter_tolerance == 1.0e-10
    assert specification.rank_threshold_multiplier == 100.0
    assert specification.start_agreement_max_abs == 1.0e-5

@pytest.mark.campaign
def test_water_factor_fit_is_rank_one_and_start_confirmed() -> None:
    result = epcsaft_regression.fit_figiel_water_solvation_factor()

    assert result.specification_id == (
        epcsaft_regression.FIGIEL_WATER_SOLVATION_FACTOR_V1.specification_id
    )
    assert tuple(start.name for start in result.starts) == ("primary", "upper")
    assert len(result.input_row_ids) == 21
    assert result.evaluated_row_ids == result.input_row_ids
    assert result.failed_row_ids == ()
    assert result.solver_converged
    assert result.numerically_converged
    assert result.physically_valid
    assert result.workflow_valid
    assert result.predictive_status == (
        "NOT_ADJUDICATED_NO_APPROVED_HELD_OUT_CUTOFF"
    )
    assert math.isfinite(result.fitted_water_solvation_factor)
    assert result.fitted_water_solvation_factor == pytest.approx(
        1.5590515389548207, rel=1.0e-12, abs=1.0e-12
    )
    assert 1.0 < result.fitted_water_solvation_factor < 2.0
    assert result.start_parameter_max_abs_delta <= 1.0e-5
    assert result.failure_reasons == ()

    for start in result.starts:
        assert start.termination == "CONVERGENCE"
        assert start.solution_usable
        assert start.final_cost <= start.initial_cost
        assert start.rank == 1
        assert start.condition_number == 1.0
        assert start.complete_jacobian_column
        assert not start.active_bound
        assert len(start.rows) == 21
        assert all(
            math.isfinite(row.modeled_gamma_pm_m)
            and row.modeled_gamma_pm_m > 0.0
            and math.isfinite(row.exact_scaled_residual_derivative)
            and row.exact_scaled_residual_derivative
            == -(
                row.modeled_gamma_pm_m / row.observed_gamma_pm_m
            )
            * row.provider_log_derivative
            for row in start.rows
        )


@pytest.mark.campaign
def test_water_factor_exact_jacobian_matches_centered_callback_values() -> None:
    specification = epcsaft_regression.FIGIEL_WATER_SOLVATION_FACTOR_V1
    model = _fixed_water_factor_model(specification)
    payload = _water_factor_native_payload(specification)
    capsule = native_sdk(model)
    trial = 1.4
    center = _native.evaluate_figiel_water_factor(capsule, payload, trial)
    differences = []
    for step in (1.0e-4, 5.0e-5):
        plus = _native.evaluate_figiel_water_factor(
            capsule, payload, trial + step
        )
        minus = _native.evaluate_figiel_water_factor(
            capsule, payload, trial - step
        )
        differences.append(
            tuple(
                (plus_row[5] - minus_row[5]) / (2.0 * step)
                for plus_row, minus_row in zip(plus, minus, strict=True)
            )
        )

    for index, row in enumerate(center):
        exact = row[7]
        tolerance = max(
            1.0e-8,
            20.0 * abs(differences[0][index] - differences[1][index]),
            2.0e-8 * abs(exact),
        )
        assert exact == pytest.approx(
            differences[1][index], rel=0.0, abs=tolerance
        )
