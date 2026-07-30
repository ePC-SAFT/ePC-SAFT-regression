from __future__ import annotations

import json
from dataclasses import replace

import pytest
from epcsaft import Mixture, Parameters

from epcsaft_regression import (
    AcquisitionClass,
    AffineParameterTransform,
    ConfirmationControls,
    CorrelationProvenance,
    FixedCompositionVleObservation,
    ObjectiveContract,
    ObservationDataset,
    PairParameterIdentity,
    ParameterFamily,
    ParameterRequest,
    RankControls,
    ReproductionClass,
    ResultContext,
    RowProvenance,
    SolverControls,
    SourceInput,
    prepare_fit,
    usability,
)

SHA = "a" * 64


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


def _model() -> Mixture:
    return Mixture(
        Parameters.from_catalog(
            "gross-2001-methane-ethane",
            components=("methane", "ethane"),
            version=1,
        )
    )


def _prepared():
    model = _model()
    request = ParameterRequest(
        family=ParameterFamily.K_IJ,
        identity=PairParameterIdentity("methane", "ethane"),
        transform=AffineParameterTransform(0.0, 0.01),
        lower_bound=-0.15,
        upper_bound=0.10,
    )
    return prepare_fit(
        model,
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


def test_dataset_records_bind_provenance_and_reject_unknown_fields() -> None:
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

    correlation = CorrelationProvenance(
        "a + b/T",
        (("a", 1.0), ("b", 2.0)),
        "Pa",
        "200 K <= T <= 300 K",
        (200.0, 250.0, 300.0),
        "sampled exactly at the declared grid",
    )
    correlated = replace(
        dataset.row_provenance[0][1],
        acquisition=AcquisitionClass.AUTHOR_CORRELATION,
        correlation=correlation,
    )
    assert correlated.correlation is correlation

    records = (dict(_records()[0], unexpected="not allowed"),)

    with pytest.raises(ValueError, match=r"row-1.*unexpected"):
        ObservationDataset.from_records(
            FixedCompositionVleObservation,
            records,
            source=_source(),
            objective=_objective(),
            row_provenance={"row-1": _dataset().row_provenance[0][1]},
        )
    with pytest.raises(ValueError, match="row IDs must be unique"):
        ObservationDataset.from_records(
            FixedCompositionVleObservation,
            _records() * 2,
            source=_source(),
            objective=_objective(),
            row_provenance={"row-1": dataset.row_provenance[0][1]},
        )


def test_prepare_and_preflight_resolve_capability_and_exact_rank(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    prepared = _prepared()
    problem = prepared.problem
    support = usability.support_view(prepared.model)
    capability = next(
        item
        for item in support
        if getattr(item, "family", None) is ParameterFamily.K_IJ
    )

    assert problem.parameters[0].capability_id == capability.capability_id
    assert problem.parameters[0].provider_parameter_fingerprint == (
        capability.parameter_fingerprint
    )
    assert problem.parameters[0].unit == "1"
    assert problem.maximum_iterations == 200
    assert problem.maximum_condition_number == 1.0e10
    assert capability.installed_ready
    assert capability.unsupported_reason == ""

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


@pytest.fixture(scope="module")
def installed_fit():
    prepared = _prepared()
    return prepared, prepared.preflight(), prepared.fit()


def test_installed_prepare_fit_and_record_contract(installed_fit) -> None:
    prepared, report, result = installed_fit
    assert report.ready
    assert report.starts[0].full_rank == result.jacobian.full_rank
    assert report.starts[0].projected_parameter_rank == (
        result.jacobian.projected_parameter_rank
    )
    assert result.workflow_valid

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
    assert record["parameters"][0]["final"] == result.parameters[0].final
    assert record["literature"]["source_printed_parameters"]["k_ij"] == "-0.01"

    with pytest.raises(ValueError, match="finite"):
        replace(result, final_cost=float("nan")).to_json_bytes(prepared=prepared)
