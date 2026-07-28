from __future__ import annotations

import math
from dataclasses import replace

from epcsaft import EPCSAFT, ParameterBundle
import pytest

from epcsaft_regression import (
    AffineParameterTransform,
    FixedCompositionVleObservation,
    ObservationPartition,
    PairParameterIdentity,
    ParameterCoordinate,
    ParameterFamily,
    RegressionProblem,
    RegressionResult,
    SourceDescriptor,
    canonical_dataset_sha256,
    fit_parameters,
    parameter_capabilities,
)
from epcsaft_regression.parameter_regression import _evaluate_parameters


def _model() -> EPCSAFT:
    parameters = ParameterBundle.from_catalog(
        "gross-2001-methane-ethane", version=1
    ).select(("methane", "ethane"))
    return EPCSAFT(parameters)


def _row(
    row_id: str = "may-2015-row-1",
    partition: ObservationPartition = ObservationPartition.TRAINING,
) -> FixedCompositionVleObservation:
    return FixedCompositionVleObservation(
        row_id=row_id,
        source_id="may-2015",
        source_locator=f"table:{row_id}",
        component_ids=("methane", "ethane"),
        temperature_k=223.15,
        pressure_pa=2.0e6,
        liquid_mole_fraction_first=0.35,
        vapor_mole_fraction_first=0.80,
        pressure_scale_pa=2.0e6,
        chemical_potential_scales=(1.0, 1.0),
        liquid_volume_origin_m3_per_mol=6.0e-5,
        liquid_volume_start_m3_per_mol=6.5e-5,
        liquid_volume_bounds_m3_per_mol=(2.0e-5, 1.0e-4),
        vapor_volume_origin_m3_per_mol=9.0e-4,
        vapor_volume_start_m3_per_mol=1.0e-3,
        vapor_volume_bounds_m3_per_mol=(1.0e-4, 1.0e-2),
        partition=partition,
    )


def _problem(
    model: EPCSAFT,
    rows: tuple[FixedCompositionVleObservation, ...] | None = None,
) -> RegressionProblem:
    (capability,) = parameter_capabilities(model)
    observations = rows or (_row(),)
    source = SourceDescriptor(
        source_id="may-2015",
        citation="Audited May 2015 methane/ethane VLE row.",
        durable_locator="local-audited-source-packet",
        source_artifact_sha256="a" * 64,
        canonical_dataset_sha256=canonical_dataset_sha256(observations),
        transformation_record="No transformation for this contract test.",
        units_and_bases="T/K, P/Pa, mole fractions.",
        use_basis="Regression derivative contract evidence.",
        residual_scale_rationale="Pressure by observed P; mu/RT dimensionless.",
    )
    parameter = ParameterCoordinate(
        family=ParameterFamily.K_IJ,
        identity=PairParameterIdentity("methane", "ethane"),
        capability_id=capability.capability_id,
        provider_parameter_fingerprint=capability.parameter_fingerprint,
        provider_topology_fingerprint=capability.topology_fingerprint,
        unit="1",
        transform=AffineParameterTransform(origin=0.0, scale=0.01),
        lower_bound=-0.15,
        upper_bound=0.10,
        starts=(0.0, -0.05, 0.05),
    )
    return RegressionProblem(
        sources=(source,),
        parameters=(parameter,),
        observations=observations,
        maximum_condition_number=1.0e10,
        maximum_iterations=50,
        function_tolerance=1.0e-12,
        gradient_tolerance=1.0e-12,
        parameter_tolerance=1.0e-12,
        confirmation_parameter_scaled_max_delta=1.0e-5,
        confirmation_cost_relative_delta=1.0e-8,
    )


def test_installed_provider_advertises_exact_neutral_binary_kij_contract() -> None:
    model = _model()

    (capability,) = parameter_capabilities(model)

    assert capability.capability_id == "neutral_binary_phase_kij_v1"
    assert capability.family is ParameterFamily.K_IJ
    assert capability.component_ids == ("methane", "ethane")
    assert capability.coordinate_kinds == ("amount", "amount", "volume", "k_ij")
    assert capability.coordinate_units == ("mol", "mol", "m3", "dimensionless")
    assert capability.parameter_fingerprint == model.parameter_fingerprint
    assert capability.topology_fingerprint.startswith("sha256:")
    assert capability.derivative_order == 2
    assert capability.maturity == "DERIVATIVE_READY"
    assert capability.authority_effect == "NONE"


def test_exact_lifted_kij_jacobian_matches_directional_residual_difference() -> None:
    model = _model()
    problem = _problem(model)
    variables = (
        0.0,
        math.log(6.5e-5 / 6.0e-5),
        math.log(1.0e-3 / 9.0e-4),
    )
    direction = (0.20, -0.10, 0.15)

    residuals, jacobian = _evaluate_parameters(problem, model, variables)
    step = 1.0e-6
    plus = tuple(
        value + step * delta
        for value, delta in zip(variables, direction, strict=True)
    )
    minus = tuple(
        value - step * delta
        for value, delta in zip(variables, direction, strict=True)
    )
    residuals_plus = _evaluate_parameters(problem, model, plus)[0]
    residuals_minus = _evaluate_parameters(problem, model, minus)[0]
    finite_difference = tuple(
        (upper - lower) / (2.0 * step)
        for upper, lower in zip(residuals_plus, residuals_minus, strict=True)
    )
    exact_product = tuple(
        math.fsum(jacobian[row * 3 + column] * direction[column] for column in range(3))
        for row in range(4)
    )

    assert len(residuals) == 4
    assert len(jacobian) == 12
    assert exact_product == pytest.approx(finite_difference, rel=2.0e-7, abs=2.0e-8)


def test_general_kij_fit_reports_rank_confirmation_and_partition_isolation() -> None:
    model = _model()
    training = _row()
    held = replace(
        _row("may-2015-held", ObservationPartition.HELD_OUT),
        pressure_pa=2.2e6,
        pressure_scale_pa=2.2e6,
    )

    result = fit_parameters(_problem(model, (training, held)), model)
    perturbed_held = replace(held, pressure_pa=8.0e6, pressure_scale_pa=8.0e6)
    repeated = fit_parameters(
        _problem(model, (training, perturbed_held)), model
    )

    assert isinstance(result, RegressionResult)
    assert result.parameter.final == repeated.parameter.final
    assert result.final_cost == repeated.final_cost
    assert result.parameter.lower_bound <= result.parameter.final <= result.parameter.upper_bound
    assert result.jacobian.full_rank == 3
    assert result.jacobian.projected_parameter_rank == 1
    assert result.confirmation_count == 2
    assert result.training_row_count == 1
    assert result.held_out_row_count == 1
    assert result.stress_row_count == 0
    assert tuple(row.partition for row in result.rows) == ("training", "held_out")
    assert result.predictive_status == "NOT_ADJUDICATED_NO_APPROVED_HELD_OUT_CUTOFF"
