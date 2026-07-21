from __future__ import annotations

from dataclasses import replace
import math

import pytest
from epcsaft import EPCSAFT, ParameterBundle, native_sdk

import epcsaft_regression
import epcsaft_regression._native as native
from epcsaft_regression.records import FIGIEL_BORN_DIAMETER_TRACER_V1
from epcsaft_regression.workflow import _born_native_payload


def _models() -> tuple[EPCSAFT, ...]:
    catalog = ParameterBundle.from_catalog("figiel-2025-reference-electrolytes", version=1)
    return tuple(
        EPCSAFT(catalog.select(target.component_order))
        for target in FIGIEL_BORN_DIAMETER_TRACER_V1.targets
    )


def test_figiel_tracer_is_one_closed_five_target_contract() -> None:
    specification = FIGIEL_BORN_DIAMETER_TRACER_V1

    assert specification.specification_id == "figiel-2025-five-ion-born-diameter-tracer-v1"
    assert specification.target_ids == (
        "figiel2025-s5-Lip-reported-average",
        "figiel2025-s5-Nap-reported-average",
        "figiel2025-s5-Kp-reported-average",
        "figiel2025-s5-Clm-reported-average",
        "figiel2025-s5-Brm-reported-average",
    )
    assert specification.ion_labels == ("Li+", "Na+", "K+", "Cl-", "Br-")
    assert specification.active_component_ids == (
        "lithium-cation",
        "sodium-cation",
        "potassium-cation",
        "chloride-anion",
        "bromide-anion",
    )
    assert specification.targets_j_per_mol == (
        -486_200.0,
        -381_100.0,
        -309_100.0,
        -314_900.0,
        -290_900.0,
    )
    assert specification.published_diameters_angstrom == (2.784, 3.445, 4.150, 4.100, 4.480)
    assert specification.temperature_k == 298.15
    assert specification.pressure_pa == 100_000.0
    assert specification.start_diameters_angstrom == (
        (3.0,) * 5,
        (2.0,) * 5,
        (5.0,) * 5,
    )
    assert specification.diameter_bounds_angstrom == (1.0, 6.0)
    assert specification.diameter_origin_angstrom == 3.0
    assert specification.diameter_scale_angstrom == 1.0
    assert specification.scaled_bounds == (-2.0, 3.0)
    assert specification.max_num_iterations == 500
    assert specification.function_tolerance == 1.0e-10
    assert specification.gradient_tolerance == 1.0e-10
    assert specification.parameter_tolerance == 1.0e-10
    assert specification.scaled_residual_max == 1.0e-8
    assert specification.confirmation_parameter_scaled_max_delta == 1.0e-5
    assert specification.observable_round_trip_j_per_mol == 50.0
    assert specification.published_diameter_reporting_half_increment_angstrom == 0.0005


def test_figiel_tracer_rejects_target_or_scope_mutation() -> None:
    with pytest.raises(ValueError, match="exact five-target contract"):
        replace(
            FIGIEL_BORN_DIAMETER_TRACER_V1,
            targets=(
                replace(
                    FIGIEL_BORN_DIAMETER_TRACER_V1.targets[0],
                    target_j_per_mol=-486_201.0,
                ),
                *FIGIEL_BORN_DIAMETER_TRACER_V1.targets[1:],
            ),
        )


def test_public_surface_adds_only_the_closed_figiel_workflow() -> None:
    assert epcsaft_regression.FIGIEL_BORN_DIAMETER_TRACER_V1 is FIGIEL_BORN_DIAMETER_TRACER_V1
    assert hasattr(epcsaft_regression, "BornDiameterFitResult")
    assert hasattr(epcsaft_regression, "fit_figiel_born_diameters")
    for excluded in (
        "ParameterRegistry",
        "ParameterOverlay",
        "fit_electrolyte_parameters",
        "persist_provider_parameters",
    ):
        assert not hasattr(epcsaft_regression, excluded)


def test_installed_provider_born_derivatives_match_step_halved_value_differences() -> None:
    specification = FIGIEL_BORN_DIAMETER_TRACER_V1
    models = _models()
    capsules = tuple(native_sdk(model) for model in models)
    payload = _born_native_payload(specification)
    trial = (2.9, 3.2, 4.3, 3.9, 4.7)
    residuals, jacobian, rows, fingerprints, compiled_identity = native.evaluate_born(
        capsules, payload, trial
    )

    assert tuple(fingerprints) == tuple(
        target.expected_provider_fingerprint for target in specification.targets
    )
    assert tuple(compiled_identity) == payload[0]
    assert len(residuals) == 5
    assert len(jacobian) == 25
    assert all(jacobian[row * 5 + column] == 0.0 for row in range(5) for column in range(5) if row != column)
    for index, (diameter, row) in enumerate(zip(trial, rows, strict=True)):
        step = 1.0e-4
        values = []
        for h in (step, step / 2.0):
            plus = list(trial)
            minus = list(trial)
            plus[index] = diameter + h
            minus[index] = diameter - h
            value_plus = native.evaluate_born(capsules, payload, tuple(plus))[2][index][0]
            value_minus = native.evaluate_born(capsules, payload, tuple(minus))[2][index][0]
            values.append((value_plus - value_minus) / (2.0 * h))
        derivative = float(row[1])
        tolerance = max(1.0e-3, 20.0 * abs(values[0] - values[1]), 2.0e-8 * abs(derivative))
        assert abs(derivative - values[1]) <= tolerance
        assert jacobian[index * 5 + index] == pytest.approx(
            derivative / abs(specification.targets_j_per_mol[index]), rel=0.0, abs=0.0
        )


def test_five_ion_born_fit_accepts_observable_recovery_and_reports_parameter_deltas() -> None:
    result = epcsaft_regression.fit_figiel_born_diameters(models=_models())

    assert result.solver_converged
    assert result.numerically_converged
    assert result.workflow_valid
    assert result.scientifically_valid
    assert result.predictive_status == "NOT_ADJUDICATED_NO_APPROVED_HELD_OUT_CUTOFF"
    assert tuple(start.name for start in result.starts) == ("primary", "lower", "upper")
    assert all(start.termination == "CONVERGENCE" for start in result.starts)
    assert all(start.solution_usable for start in result.starts)
    assert all(start.rank == 5 for start in result.starts)
    assert all(math.isfinite(start.condition_number) for start in result.starts)
    assert all(not parameter.active_bound for parameter in result.parameters)
    assert tuple(parameter.final_diameter_angstrom for parameter in result.parameters) == pytest.approx(
        (2.7888130173797934, 3.4524616464076425, 4.147266741279482, 4.101505615791675, 4.476998527506598),
        rel=2.0e-12,
        abs=2.0e-12,
    )
    assert all(
        abs(row.raw_error_j_per_mol)
        <= FIGIEL_BORN_DIAMETER_TRACER_V1.observable_round_trip_j_per_mol
        for row in result.observations
    )
    assert result.confirmation_parameter_scaled_max_deltas[0] <= 1.0e-5
    assert result.confirmation_parameter_scaled_max_deltas[1] <= 1.0e-5
    # The published Table 3 values are rounded comparison anchors, not residual
    # targets or a second acceptance gate.
    assert max(abs(parameter.published_delta_angstrom) for parameter in result.parameters) > 0.0005
    assert result.failure_reasons == ()


def test_born_workflow_rejects_wrong_model_order_before_ceres() -> None:
    models = _models()
    with pytest.raises(ValueError, match="component order"):
        epcsaft_regression.fit_figiel_born_diameters(models=(models[1], models[0], *models[2:]))
