from __future__ import annotations

from dataclasses import replace

import pytest

from epcsaft_regression import (
    AffineParameterTransform,
    ComponentParameterIdentity,
    EvaluatorContract,
    ObservationPartition,
    ParameterCoordinate,
    ParameterFamily,
    PositiveEvaluatorProblem,
    PositiveObservationTransform,
    PositiveScalarObservation,
    SourceDescriptor,
    canonical_positive_dataset_sha256,
)


def _rows() -> tuple[PositiveScalarObservation, ...]:
    return (
        PositiveScalarObservation(
            row_id="pressure-1",
            state_id="state-1",
            state_schema_id="fixed-tp-liquid-v1",
            source_id="doi:example",
            source_locator="table-1:pressure-1",
            primitive_id="neutral_component_fugacity_pa;co2",
            primitive_unit="Pa",
            transform=PositiveObservationTransform.NATURAL_LOG,
            reference_id="provider-helmholtz-coordinate-basis",
            reference_fingerprint="sha256:" + "1" * 64,
            observed_value=574.0,
            residual_scale=2.302585092994046,
            residual_scale_unit="1",
            partition=ObservationPartition.TRAINING,
        ),
        PositiveScalarObservation(
            row_id="speciation-1",
            state_id="state-1",
            state_schema_id="fixed-tp-liquid-v1",
            source_id="doi:example",
            source_locator="table-2:speciation-1",
            primitive_id="species_mole_fraction;carbamate",
            primitive_unit="dimensionless",
            transform=PositiveObservationTransform.IDENTITY,
            reference_id="source-standard-state-v1",
            reference_fingerprint="sha256:" + "2" * 64,
            observed_value=0.0502,
            residual_scale=0.01,
            residual_scale_unit="dimensionless",
            partition=ObservationPartition.TRAINING,
        ),
    )


def _contract() -> EvaluatorContract:
    return EvaluatorContract(
        evaluator_identity="owner.homogeneous-liquid-observation.v1",
        capability_id="homogeneous-liquid-positive-scalars-v1",
        capability_fingerprint="sha256:" + "3" * 64,
        provider_artifact_identity=(
            "epcsaft==0.2.0.dev0;RECORD=sha256:"
            + "4" * 64
            + ";HEADER=sha256:"
            + "5" * 64
        ),
        owner_artifact_identity=(
            "owner==0.2.0.dev0;RECORD=sha256:" + "6" * 64 + ";HEADER=sha256:" + "7" * 64
        ),
        contract_fingerprint="sha256:" + "8" * 64,
        artifact_identity="sha256:" + "c" * 64,
    )


def _problem() -> PositiveEvaluatorProblem:
    rows = _rows()
    source = SourceDescriptor(
        source_id="doi:example",
        citation="Example et al. (2026)",
        durable_locator="https://doi.org/10.0000/example",
        source_artifact_sha256="d" * 64,
        canonical_dataset_sha256=canonical_positive_dataset_sha256(rows),
        transformation_record="No transformation beyond declared residual transforms.",
        units_and_bases="Pressure in Pa; speciation on mole-fraction basis.",
        use_basis="Source-bound transport contract test.",
        residual_scale_rationale="Declared engineering scales; not uncertainty.",
    )
    parameters = tuple(
        ParameterCoordinate(
            family=ParameterFamily.SEGMENT_DIAMETER,
            identity=ComponentParameterIdentity(component_id),
            capability_id="reacting_phase_active_parameter_v1",
            provider_parameter_fingerprint="sha256:" + "a" * 64,
            provider_topology_fingerprint="sha256:" + "b" * 64,
            unit="angstrom",
            transform=AffineParameterTransform(origin=3.0, scale=1.0),
            lower_bound=2.0,
            upper_bound=5.8,
        )
        for component_id in ("cation", "carbamate")
    )
    return PositiveEvaluatorProblem(
        sources=(source,),
        parameters=parameters,
        parameter_ids=(
            "segment_diameter;component;cation",
            "segment_diameter;component;carbamate",
        ),
        start_vectors=((3.4, 3.5), (2.5, 4.5)),
        observations=rows,
        evaluator=_contract(),
        maximum_condition_number=1.0e6,
        maximum_iterations=100,
        maximum_solver_time_seconds=30.0,
        function_tolerance=1.0e-12,
        gradient_tolerance=1.0e-12,
        parameter_tolerance=1.0e-12,
        confirmation_parameter_scaled_max_delta=1.0e-6,
        confirmation_cost_relative_delta=1.0e-8,
    )


def test_positive_evaluator_problem_binds_sources_rows_parameters_and_transforms() -> (
    None
):
    problem = _problem()

    assert problem.training_observations == problem.observations
    assert problem.solver_start_vectors[0] == pytest.approx((0.4, 0.5))
    assert problem.solver_start_vectors[1] == pytest.approx((-0.5, 1.5))
    assert problem.parameter_ids == (
        "segment_diameter;component;cation",
        "segment_diameter;component;carbamate",
    )
    assert problem.observations[0].transform is PositiveObservationTransform.NATURAL_LOG
    assert problem.observations[1].transform is PositiveObservationTransform.IDENTITY


@pytest.mark.parametrize(
    ("index", "changes"),
    (
        (0, {"observed_value": 0.0}),
        (0, {"residual_scale_unit": "Pa"}),
        (1, {"residual_scale_unit": "1"}),
    ),
)
def test_positive_observation_rejects_nonpositive_or_dimensionally_invalid_rows(
    index: int,
    changes: dict[str, object],
) -> None:
    with pytest.raises(ValueError):
        replace(_rows()[index], **changes)


def test_positive_evaluator_problem_rejects_parameter_or_fingerprint_mismatch() -> None:
    with pytest.raises(ValueError, match="parameter_ids"):
        replace(_problem(), parameter_ids=("segment_diameter;component;cation",))

    with pytest.raises(ValueError, match="Provider parameter fingerprint"):
        replace(
            _problem(),
            parameters=(
                _problem().parameters[0],
                replace(
                    _problem().parameters[1],
                    provider_parameter_fingerprint="sha256:" + "e" * 64,
                ),
            ),
        )
