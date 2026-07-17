from __future__ import annotations

import math

import pytest
from epcsaft import EPCSAFT, ParameterBundle, native_sdk

import epcsaft_regression._native as native
from epcsaft_regression.records import METHANE_FIT_SPECIFICATION_V1, load_methane_dataset
from epcsaft_regression.workflow import _native_payload, fit_methane_saturation


def _capsule() -> object:
    parameters = ParameterBundle.from_catalog(
        "gross-2001-methane-ethane", version=1
    ).select(("methane",))
    return native_sdk(EPCSAFT(parameters))


def _payload() -> tuple[object, ...]:
    dataset = load_methane_dataset()
    specification = METHANE_FIT_SPECIFICATION_V1
    model = _model()
    return _native_payload(dataset, specification, model.parameter_fingerprint)


def _model() -> EPCSAFT:
    parameters = ParameterBundle.from_catalog(
        "gross-2001-methane-ethane", version=1
    ).select(("methane",))
    return EPCSAFT(parameters)


def test_start_residuals_match_provider_direct_130_k_anchor() -> None:
    residuals, jacobian, diagnostics, fingerprint = native.evaluate(
        _capsule(), _payload(), (0.0,) * 11
    )

    assert residuals[4:8] == pytest.approx(
        (
            -30.793189605316272,
            -0.043171211223460487,
            -0.64887647749645083,
            0.0,
        ),
        rel=2.0e-13,
        abs=2.0e-13,
    )
    assert len(residuals) == 16
    assert len(jacobian) == 16 * 11
    assert len(diagnostics) == 4
    assert fingerprint.startswith("sha256:")


def test_exact_jacobian_matches_independent_directional_residual_difference() -> None:
    capsule = _capsule()
    payload = _payload()
    variables = (0.0,) * 11
    direction = (
        0.2,
        -0.1,
        0.05,
        0.01,
        -0.02,
        -0.015,
        0.012,
        0.008,
        -0.01,
        -0.006,
        0.014,
    )
    residuals, jacobian, _, _ = native.evaluate(capsule, payload, variables)
    step = 1.0e-6
    plus = tuple(value + step * delta for value, delta in zip(variables, direction, strict=True))
    minus = tuple(value - step * delta for value, delta in zip(variables, direction, strict=True))
    residuals_plus = native.evaluate(capsule, payload, plus)[0]
    residuals_minus = native.evaluate(capsule, payload, minus)[0]
    finite_difference = tuple(
        (right - left) / (2.0 * step)
        for right, left in zip(residuals_plus, residuals_minus, strict=True)
    )
    product = tuple(
        math.fsum(jacobian[row * 11 + column] * direction[column] for column in range(11))
        for row in range(16)
    )

    assert residuals
    assert product == pytest.approx(finite_difference, rel=2.0e-6, abs=2.0e-7)


def test_public_workflow_returns_strict_candidate_diagnostics() -> None:
    parameters = ParameterBundle.from_catalog(
        "gross-2001-methane-ethane", version=1
    ).select(("methane",))
    model = EPCSAFT(parameters)
    dataset = load_methane_dataset()

    result = fit_methane_saturation(
        model=model,
        dataset=dataset,
        specification=METHANE_FIT_SPECIFICATION_V1,
    )

    assert result.solver_converged, result.failure_reasons
    assert result.numerically_converged, result.failure_reasons
    assert result.physically_valid, result.failure_reasons
    assert result.termination == "CONVERGENCE"
    assert result.solution_usable
    assert math.isfinite(result.initial_cost)
    assert math.isfinite(result.final_cost)
    assert result.final_cost < result.initial_cost
    assert result.provider_fingerprint == model.parameter_fingerprint
    assert result.provider_fingerprint in result.compiled_problem_identity
    assert (
        "f92f79c8d6f614660e5c201b7061c9b02b5cd1a25a4ed8c8fee0b59adaabf2bf"
        in result.compiled_problem_identity
    )
    assert len(result.parameters) == 3
    assert any(abs(item.movement) > 1.0e-8 for item in result.parameters)
    assert all(item.lower_bound <= item.final <= item.upper_bound for item in result.parameters)
    assert not any(item.active_bound for item in result.parameters)
    assert result.jacobian.complete_columns
    assert result.jacobian.full_rank == 11
    assert len(result.jacobian.full_singular_values) == 11
    assert result.jacobian.parameter_rank == 3
    assert len(result.jacobian.parameter_singular_values) == 3
    assert math.isfinite(result.jacobian.full_condition_number)
    assert math.isfinite(result.jacobian.parameter_condition_number)
    assert len(result.training_rows) == 4
    assert all(len(row.raw_residuals) == 4 for row in result.training_rows)
    assert all(len(row.scaled_residuals) == 4 for row in result.training_rows)
    assert all(row.liquid_volume_m3 < row.vapor_volume_m3 for row in result.training_rows)
    assert all(row.liquid_stability_slope > 0.0 for row in result.training_rows)
    assert all(row.vapor_stability_slope > 0.0 for row in result.training_rows)
    assert len(result.reporting_rows) == 9
    assert tuple(row.temperature_k for row in result.reporting_rows) == tuple(
        float(value) for value in range(100, 181, 10)
    )
    assert sum(row.training for row in result.reporting_rows) == 4
    assert all(row.physically_valid for row in result.reporting_rows)
    assert all(row.predicted_pressure_pa > 0.0 for row in result.reporting_rows)
    assert all(row.predicted_liquid_density_kg_m3 > 0.0 for row in result.reporting_rows)
    assert result.confirmation_parameter_scaled_max_delta <= 1.0e-5
    assert result.confirmation_cost_relative_delta <= 1.0e-8
    assert result.failure_reasons == ()
