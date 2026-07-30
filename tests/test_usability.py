from __future__ import annotations

from dataclasses import replace
import json
import math

import pytest

import epcsaft_regression.usability as usability
from epcsaft import Mixture, Parameters
from epcsaft_regression import (
    AcquisitionClass,
    AffineParameterTransform,
    ComponentParameterIdentity,
    ConfirmationControls,
    FixedCompositionVleObservation,
    ObjectiveContract,
    ObservationDataset,
    ObservationPartition,
    ParameterCapability,
    ParameterFamily,
    ParameterRequest,
    PairParameterIdentity,
    RankControls,
    ReproductionClass,
    ResultContext,
    RowProvenance,
    SolverControls,
    SourceInput,
    prepare_fit,
)
from epcsaft_regression.parameter_regression import (
    FittedParameterDiagnostic,
    GeneralJacobianDiagnostics,
    GeneralRowDiagnostic,
    RegressionResult,
)


SHA = "a" * 64
PREFIXED_SHA = "sha256:" + "b" * 64


def _source() -> SourceInput:
    return SourceInput(
        source_id="doi:example",
        citation="Example et al. (2026)",
        durable_locator="https://doi.org/10.0000/example",
        source_artifact_sha256=SHA,
        transformation_record="bar to Pa; no exclusions",
        units_and_bases="T/K, P/Pa, mole fractions",
        use_basis="source-backed test data",
        residual_scale_rationale="pressure by observed pressure; mu/RT",
    )


def _objective() -> ObjectiveContract:
    return ObjectiveContract(
        residual_family="fixed_composition_vle",
        interpretation="native_scaled_least_squares",
        row_weighting="scales embedded in each validated row",
        covariance_interpretation="independent rows; no covariance supplied",
        loss="squared",
        loss_parameters=(),
        failed_row_treatment="fail",
    )


def _records() -> tuple[dict[str, object], ...]:
    return (
        {
            "row_id": "row-1",
            "source_id": "doi:example",
            "source_locator": "table-1:row-1",
            "component_ids": ["methane", "ethane"],
            "temperature_k": 203.22,
            "pressure_pa": 2_124_000.0,
            "liquid_mole_fraction_first": 0.3653,
            "vapor_mole_fraction_first": 0.8667,
            "pressure_scale_pa": 2_124_000.0,
            "chemical_potential_scales": [1.0, 1.0],
            "liquid_volume_origin_m3_per_mol": 6.0e-5,
            "liquid_volume_start_m3_per_mol": 6.5e-5,
            "liquid_volume_bounds_m3_per_mol": [2.0e-5, 1.0e-4],
            "vapor_volume_origin_m3_per_mol": 9.0e-4,
            "vapor_volume_start_m3_per_mol": 1.0e-3,
            "vapor_volume_bounds_m3_per_mol": [1.0e-4, 1.0e-2],
            "partition": "training",
        },
    )


def _dataset() -> ObservationDataset:
    return ObservationDataset.from_records(
        FixedCompositionVleObservation,
        _records(),
        source=_source(),
        objective=_objective(),
        row_provenance={
            "row-1": RowProvenance(
                acquisition=AcquisitionClass.DIRECT_MEASUREMENT,
                duplicate_decision="unique source row",
                exclusion_decision="included",
                critical_region_decision="outside critical region",
                censoring_decision="not censored",
                outlier_decision="retained; no outlier rule applied",
            )
        },
    )


def _capability() -> ParameterCapability:
    return ParameterCapability(
        capability_id="neutral_binary_phase_kij_v1",
        family=ParameterFamily.K_IJ,
        component_ids=("methane", "ethane"),
        coordinate_kinds=("amount", "amount", "volume", "k_ij"),
        coordinate_units=("mol", "mol", "m3", "dimensionless"),
        parameter_fingerprint=PREFIXED_SHA,
        topology_fingerprint="sha256:" + "c" * 64,
        derivative_order=2,
        maturity="DERIVATIVE_READY",
        authority_effect="none",
        temperature_min_k=100.0,
        temperature_max_k=400.0,
        identity_shape="unordered_component_pair",
        observation_contract="fixed_composition_helmholtz_phase",
        model_domain="neutral_nonassociating_binary",
        tensor_layout="row_major",
        state_coordinate_count=3,
        active_parameter_count=1,
        helmholtz_basis_id="molar_helmholtz",
        unsupported_status="",
        domain_status="ready",
        active_component_ids=("ethane", "methane"),
    )


def _prepared(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(usability, "parameter_capabilities", lambda model: (_capability(),))
    request = ParameterRequest(
        family=ParameterFamily.K_IJ,
        identity=PairParameterIdentity("methane", "ethane"),
        transform=AffineParameterTransform(0.0, 0.01),
        lower_bound=-0.15,
        upper_bound=0.10,
    )
    return prepare_fit(
        object(),
        datasets=(_dataset(),),
        parameters=(request,),
        parameter_slot_indices=(0,),
        start_vectors=((0.0,), (-0.05,)),
        solver=SolverControls(
            maximum_iterations=200,
            maximum_solver_time_seconds=30.0,
            function_tolerance=1.0e-12,
            gradient_tolerance=1.0e-12,
            parameter_tolerance=1.0e-12,
        ),
        rank=RankControls(maximum_condition_number=1.0e10),
        confirmation=ConfirmationControls(
            parameter_scaled_max_delta=1.0e-5,
            cost_relative_delta=1.0e-8,
        ),
    )


def test_records_construct_existing_observations_and_bind_provenance() -> None:
    dataset = _dataset()

    assert isinstance(dataset.observations[0], FixedCompositionVleObservation)
    assert dataset.observations[0].component_ids == ("methane", "ethane")
    assert dataset.row_provenance[0][1].acquisition is AcquisitionClass.DIRECT_MEASUREMENT
    assert len(dataset.provenance_sha256) == 64

    changed = replace(
        dataset.row_provenance[0][1],
        exclusion_decision="included after source-author clarification",
    )
    repeated = ObservationDataset.from_records(
        FixedCompositionVleObservation,
        _records(),
        source=_source(),
        objective=_objective(),
        row_provenance={"row-1": changed},
    )
    assert repeated.provenance_sha256 != dataset.provenance_sha256


def test_record_constructor_fails_with_row_specific_schema_reason() -> None:
    records = (dict(_records()[0], unexpected="not allowed"),)

    with pytest.raises(ValueError, match=r"row-1.*unexpected"):
        ObservationDataset.from_records(
            FixedCompositionVleObservation,
            records,
            source=_source(),
            objective=_objective(),
            row_provenance={"row-1": _dataset().row_provenance[0][1]},
        )


def test_prepare_resolves_capability_identity_and_controls(monkeypatch: pytest.MonkeyPatch) -> None:
    prepared = _prepared(monkeypatch)
    problem = prepared.problem

    assert problem.parameters[0].capability_id == _capability().capability_id
    assert problem.parameters[0].provider_parameter_fingerprint == PREFIXED_SHA
    assert problem.parameters[0].unit == "1"
    assert problem.maximum_iterations == 200
    assert problem.maximum_condition_number == 1.0e10


def test_preflight_evaluates_starts_without_solving_and_reports_rank(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    prepared = _prepared(monkeypatch)
    calls: list[tuple[float, ...]] = []

    def evaluate(problem: object, model: object, variables: tuple[float, ...]):
        calls.append(variables)
        # R=4, N=1, Q=2; all three columns are independent.
        return (
            (0.0, 0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0),
        )

    monkeypatch.setattr(usability, "_evaluate_parameters", evaluate)
    report = prepared.preflight()

    assert len(calls) == 2
    assert (report.fitted_parameter_count, report.lifted_variable_count) == (1, 2)
    assert report.residual_count == 4
    assert report.starts[0].full_rank == 3
    assert report.starts[0].projected_parameter_rank == 1
    assert report.ready
    assert report.reasons == ()


def test_installed_public_prepare_preflight_and_fit_path() -> None:
    model = Mixture(
        Parameters.from_catalog(
            "gross-2001-methane-ethane",
            components=("methane", "ethane"),
            version=1,
        )
    )
    request = ParameterRequest(
        ParameterFamily.K_IJ,
        PairParameterIdentity("methane", "ethane"),
        AffineParameterTransform(0.0, 0.01),
        -0.15,
        0.10,
    )
    prepared = prepare_fit(
        model,
        datasets=(_dataset(),),
        parameters=(request,),
        parameter_slot_indices=(0,),
        start_vectors=((0.0,), (-0.05,)),
        solver=SolverControls(50, 30.0, 1e-12, 1e-12, 1e-12),
        rank=RankControls(1e10),
        confirmation=ConfirmationControls(1e-5, 1e-8),
    )

    report = prepared.preflight()
    result = prepared.fit()

    assert report.ready
    assert report.starts[0].full_rank == result.jacobian.full_rank
    assert report.starts[0].projected_parameter_rank == (
        result.jacobian.projected_parameter_rank
    )
    assert result.workflow_valid


def _result(prepared: object) -> RegressionResult:
    problem = prepared.problem
    capability = _capability()
    return RegressionResult(
        problem=problem,
        capabilities=(capability,),
        provider_parameter_fingerprint=capability.parameter_fingerprint,
        provider_topology_fingerprint=capability.topology_fingerprint,
        solver_converged=True,
        numerically_converged=True,
        workflow_valid=True,
        physical_status="NOT_ADJUDICATED",
        scientific_status="NOT_ADJUDICATED",
        predictive_status="NOT_ADJUDICATED",
        termination="CONVERGENCE",
        solution_usable=True,
        initial_cost=1.0,
        final_cost=0.0,
        iterations=3,
        residual_evaluation_count=4,
        jacobian_evaluation_count=2,
        parameters=(
            FittedParameterDiagnostic(
                family=ParameterFamily.K_IJ,
                component_ids=("ethane", "methane"),
                unit="1",
                transform_origin=0.0,
                transform_scale=0.01,
                start=0.0,
                final=-0.008,
                movement=-0.008,
                lower_bound=-0.15,
                upper_bound=0.10,
                active_bound_distance=0.092,
                active_bound=None,
            ),
        ),
        jacobian=GeneralJacobianDiagnostics(
            residual_count=4,
            variable_count=3,
            full_singular_values=(2.0, 1.0, 0.5),
            full_rank=3,
            full_condition_number=4.0,
            projected_parameter_singular_values=(0.5,),
            projected_parameter_rank=1,
            projected_parameter_condition_number=1.0,
        ),
        rows=(
            GeneralRowDiagnostic(
                row_id="row-1",
                partition="training",
                liquid_volume_m3_per_mol=6.0e-5,
                vapor_volume_m3_per_mol=9.0e-4,
                scaled_residuals=(0.0, 0.0, 0.0, 0.0),
                observed_pressure_pa=2.0e6,
                liquid_model_pressure_pa=2.0e6,
                vapor_model_pressure_pa=2.0e6,
                chemical_potential_differences_over_rt=(0.0, 0.0),
                derivative_status="EXACT_PROVIDER_HESSIAN",
                status="evaluated",
                evaluated=True,
                failure_reason="",
            ),
        ),
        confirmation_count=1,
        confirmation_parameter_scaled_max_delta=0.0,
        confirmation_cost_relative_max_delta=0.0,
        confirmations_usable=True,
        training_row_count=1,
        held_out_row_count=0,
        stress_row_count=0,
        evaluated_row_count=1,
        skipped_row_count=0,
        failed_row_count=0,
        failure_reasons=(),
    )


def test_result_record_is_canonical_and_keeps_literature_identity(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    prepared = _prepared(monkeypatch)
    result = _result(prepared)
    context = ResultContext(
        reproduction_class=ReproductionClass.MODERN_REFIT,
        model_identity={"formulation": "PC-SAFT", "association_scheme": "none"},
        source_printed_parameters={"k_ij": "-0.01"},
        practical_identifiability_artifact=None,
        uncertainty_artifact=None,
    )

    first = result.to_json_bytes(prepared=prepared, context=context)
    second = result.to_json_bytes(prepared=prepared, context=context)
    record = json.loads(first)

    assert first == second
    assert record["schema_id"] == "epcsaft-regression-result"
    assert record["schema_version"] == 1
    assert record["literature"]["reproduction_class"] == "MODERN_REFIT"
    assert record["problem"]["datasets"][0]["rows"][0]["acquisition"] == (
        "direct_measurement"
    )
    assert record["parameters"][0]["final"] == -0.008
    assert record["literature"]["source_printed_parameters"]["k_ij"] == "-0.01"


def test_result_record_rejects_nonfinite_values(monkeypatch: pytest.MonkeyPatch) -> None:
    prepared = _prepared(monkeypatch)
    result = replace(_result(prepared), final_cost=math.nan)

    with pytest.raises(ValueError, match="finite"):
        result.to_json_bytes(prepared=prepared)
