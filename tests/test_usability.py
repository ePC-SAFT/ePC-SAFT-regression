from __future__ import annotations

import json
from dataclasses import replace
from typing import get_type_hints

import pytest
from epcsaft import Mixture, Parameters
from parameter_cases import neutral_parameters

from epcsaft_regression import (
    AcquisitionClass,
    AffineParameterTransform,
    ArtifactReference,
    ConfirmationControls,
    CorrelationProvenance,
    FixedCompositionVleObservation,
    LiteratureModelIdentity,
    ObjectiveContract,
    ObservationDataset,
    PairParameterIdentity,
    ParameterFamily,
    ParameterRequest,
    PreparedFitEvaluation,
    PureVaporPressureObservation,
    RankControls,
    RegressionResult,
    ReproductionClass,
    ResultContext,
    RowProvenance,
    SolverControls,
    SourceInput,
    canonical_dataset_sha256,
    fit_parameters,
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
        row_weighting="observation_residual_scales",
        covariance_interpretation="independent_no_covariance",
        loss="squared",
        loss_parameters=(),
        failed_row_treatment="fail_fit",
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
        Parameters.from_dictionary(neutral_parameters(("methane", "ethane")))
    )


def _prepared(dataset: ObservationDataset | None = None):
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
        datasets=(dataset or _dataset(),),
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
    assert (
        dataset.row_provenance[0][1].acquisition is AcquisitionClass.DIRECT_MEASUREMENT
    )
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
        ("temperature_k",),
        ((203.22,),),
        "sampled exactly at the declared grid",
    )
    correlated = replace(
        dataset.row_provenance[0][1],
        acquisition=AcquisitionClass.AUTHOR_CORRELATION,
        correlation=correlation,
    )
    assert correlated.correlation is correlation
    correlated_dataset = ObservationDataset.from_records(
        FixedCompositionVleObservation,
        _records(),
        source=_source(),
        objective=_objective(),
        row_provenance={"row-1": correlated},
    )
    assert correlated_dataset.provenance_sha256 != dataset.provenance_sha256
    resampled = ObservationDataset.from_records(
        FixedCompositionVleObservation,
        (dict(_records()[0], temperature_k=204.0),),
        source=_source(),
        objective=_objective(),
        row_provenance={
            "row-1": replace(
                correlated,
                correlation=replace(correlation, sampling_grid=((204.0,),)),
            )
        },
    )
    assert canonical_dataset_sha256(resampled.observations) != (
        canonical_dataset_sha256(correlated_dataset.observations)
    )
    multidimensional = replace(
        correlated,
        correlation=replace(
            correlation,
            sampling_fields=("temperature_k", "pressure_pa"),
            sampling_grid=((203.22, 2_124_000.0),),
        ),
    )
    ObservationDataset.from_records(
        FixedCompositionVleObservation,
        _records(),
        source=_source(),
        objective=_objective(),
        row_provenance={"row-1": multidimensional},
    )
    with pytest.raises(ValueError, match="sampling_grid"):
        ObservationDataset.from_records(
            FixedCompositionVleObservation,
            _records(),
            source=_source(),
            objective=_objective(),
            row_provenance={
                "row-1": replace(
                    correlated,
                    correlation=replace(correlation, sampling_grid=((200.0,),)),
                )
            },
        )
    with pytest.raises(TypeError, match="RowProvenance"):
        ObservationDataset.from_records(
            FixedCompositionVleObservation,
            _records(),
            source=_source(),
            objective=_objective(),
            row_provenance={"row-1": "not provenance"},  # type: ignore[dict-item]
        )
    with pytest.raises(TypeError, match="CorrelationProvenance"):
        replace(
            dataset.row_provenance[0][1],
            acquisition=AcquisitionClass.AUTHOR_CORRELATION,
            correlation="not correlation",  # type: ignore[arg-type]
        )
    with pytest.raises(ValueError, match="residual_family"):
        ObservationDataset.from_records(
            FixedCompositionVleObservation,
            _records(),
            source=_source(),
            objective=replace(_objective(), residual_family="pure_density"),
            row_provenance={"row-1": dataset.row_provenance[0][1]},
        )

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
    missing = dict(_records()[0])
    missing.pop("pressure_pa")
    with pytest.raises(ValueError, match=r"row-1.*missing pressure_pa"):
        ObservationDataset.from_records(
            FixedCompositionVleObservation,
            (missing,),
            source=_source(),
            objective=_objective(),
            row_provenance={"row-1": dataset.row_provenance[0][1]},
        )
    with pytest.raises(ValueError, match=r"row-1.*finite"):
        ObservationDataset.from_records(
            FixedCompositionVleObservation,
            (dict(_records()[0], temperature_k=float("nan")),),
            source=_source(),
            objective=_objective(),
            row_provenance={"row-1": dataset.row_provenance[0][1]},
        )
    with pytest.raises(ValueError, match=r"row-1.*invalid partition"):
        ObservationDataset.from_records(
            FixedCompositionVleObservation,
            (dict(_records()[0], partition="reserved"),),
            source=_source(),
            objective=_objective(),
            row_provenance={"row-1": dataset.row_provenance[0][1]},
        )
    reordered = ObservationDataset.from_records(
        FixedCompositionVleObservation,
        (dict(_records()[0], component_ids=["ethane", "methane"]),),
        source=_source(),
        objective=_objective(),
        row_provenance={"row-1": dataset.row_provenance[0][1]},
    )
    reordered_report = _prepared(reordered).preflight()
    assert not reordered_report.ready
    assert "component order" in reordered_report.reasons[0]

    vapor = ObservationDataset.from_records(
        PureVaporPressureObservation,
        (
            {
                "row_id": "vapor-1",
                "source_id": "doi:example",
                "source_locator": "table-2:row-1",
                "component_id": "methane",
                "temperature_k": 150.0,
                "pressure_pa": 1.0e6,
                "pressure_scale_pa": 1.0e6,
                "chemical_potential_scale": 1.0,
                "liquid_volume_origin_m3_per_mol": 4.0e-5,
                "liquid_volume_start_m3_per_mol": 4.0e-5,
                "liquid_volume_bounds_m3_per_mol": (2.0e-5, 8.0e-5),
                "vapor_volume_origin_m3_per_mol": 1.0e-3,
                "vapor_volume_start_m3_per_mol": 1.0e-3,
                "vapor_volume_bounds_m3_per_mol": (1.0e-4, 1.0e-2),
                "partition": "training",
            },
        ),
        source=_source(),
        objective=replace(_objective(), residual_family="pure_vapor_pressure"),
        row_provenance={"vapor-1": dataset.row_provenance[0][1]},
    )
    assert isinstance(vapor.observations[0], PureVaporPressureObservation)


@pytest.mark.parametrize(
    "acquisition",
    (
        AcquisitionClass.DIRECT_MEASUREMENT,
        AcquisitionClass.DIGITIZED_VALUE,
        AcquisitionClass.DATABASE_RECORD,
        AcquisitionClass.AUTHOR_CORRELATION,
        AcquisitionClass.RECONSTRUCTED_CORRELATION,
    ),
)
def test_every_acquisition_class_has_a_valid_provenance_contract(
    acquisition: AcquisitionClass,
) -> None:
    correlation = (
        CorrelationProvenance(
            "a + b/T",
            (("a", 1.0), ("b", 2.0)),
            "Pa",
            "200 K <= T <= 300 K",
            ("temperature_k",),
            ((203.22,),),
            "sampled exactly at the declared grid",
        )
        if acquisition
        in (
            AcquisitionClass.AUTHOR_CORRELATION,
            AcquisitionClass.RECONSTRUCTED_CORRELATION,
        )
        else None
    )

    provenance = RowProvenance(
        acquisition=acquisition,
        duplicate_decision="unique source row",
        exclusion_decision="included",
        critical_region_decision="outside critical region",
        censoring_decision="not censored",
        outlier_decision="retained; no outlier rule applied",
        correlation=correlation,
    )

    dataset = ObservationDataset.from_records(
        FixedCompositionVleObservation,
        _records(),
        source=_source(),
        objective=_objective(),
        row_provenance={"row-1": provenance},
    )
    assert dataset.row_provenance[0][1].acquisition is acquisition


@pytest.mark.parametrize(
    ("change", "match"),
    (
        ({"interpretation": "opaque"}, "interpretation"),
        ({"row_weighting": "opaque"}, "row_weighting"),
        ({"covariance_interpretation": "opaque"}, "covariance_interpretation"),
        ({"loss": "huber"}, "squared loss"),
        ({"loss_parameters": (("delta", 1.0),)}, "no loss parameters"),
        ({"failed_row_treatment": "skip"}, "failed_row_treatment"),
    ),
)
def test_objective_contract_fails_closed(change: dict[str, object], match: str) -> None:
    with pytest.raises(ValueError, match=match):
        replace(_objective(), **change)


def test_prepare_and_preflight_resolve_capability_and_exact_rank(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    prepared = _prepared()
    problem = prepared.problem
    support = usability.parameter_capabilities(prepared.model)
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
    assert report.starts[0].fitted_bounds_status == ("interior",)
    assert report.starts[0].lifted_bounds_status == ("interior", "interior")
    assert report.ready
    assert report.reasons == ()


def test_native_jacobian_diagnostics_reject_negative_dimensions() -> None:
    with pytest.raises(ValueError, match="must be positive"):
        usability._native.diagnose_jacobian(-1, 1, 2, ())
    with pytest.raises(ValueError, match="must be finite"):
        usability._native.diagnose_jacobian(1, 0, 1, (float("nan"),))


@pytest.mark.parametrize(
    ("error", "reason"),
    (
        (
            ValueError("installed capability derivative changed"),
            "unavailable_derivatives:",
        ),
        (RuntimeError("density iteration left the EOS domain"), "eos_domain_failure:"),
        (
            RuntimeError("observation payload has wrong field count"),
            "malformed_rows:",
        ),
        (TypeError("row payload shape changed"), "malformed_rows:"),
        (ValueError("component order changed"), "incompatible_contract:"),
        (ArithmeticError("nonfinite exact Jacobian"), "nonfinite_evaluation:"),
    ),
)
def test_preflight_preserves_failure_classes(
    monkeypatch: pytest.MonkeyPatch,
    error: Exception,
    reason: str,
) -> None:
    def fail(*_args: object) -> None:
        raise error

    monkeypatch.setattr(usability, "_evaluate_parameters", fail)
    report = _prepared().preflight()

    assert not report.ready
    assert all(start.failure_reason.startswith(reason) for start in report.starts)
    with pytest.raises(RuntimeError, match=f"^{reason}"):
        _prepared().evaluate((0.0,))


@pytest.mark.parametrize(
    ("residuals", "jacobian", "reason"),
    (
        ((0.0,) * 4, (0.0,) * 11, "unavailable_derivatives:"),
        ((0.0, 0.0, 0.0, float("nan")), (0.0,) * 12, "nonfinite_evaluation:"),
    ),
)
def test_public_evaluation_rejects_incomplete_or_nonfinite_payload(
    monkeypatch: pytest.MonkeyPatch,
    residuals: tuple[float, ...],
    jacobian: tuple[float, ...],
    reason: str,
) -> None:
    monkeypatch.setattr(
        usability,
        "_evaluate_parameters",
        lambda *_: (residuals, jacobian),
    )

    with pytest.raises(RuntimeError, match=f"^{reason}"):
        _prepared().evaluate((0.0,))


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
        model_identity=LiteratureModelIdentity(
            formulation="PC-SAFT",
            phase_reference_convention="homogeneous Helmholtz phases",
            association_scheme="none",
            site_multiplicities="not applicable",
            mixing_rule="one-fluid",
            combining_rule="Lorentz-Berthelot with declared k_ij",
            mixture_parameter_treatment="k_ij fitted; l_ij fixed to zero",
        ),
        reconstruction_decisions=(
            "used the current exact-EOS scaled least-squares objective",
        ),
        source_printed_parameters=(("k_ij", "-0.01"),),
        profile_artifact=ArtifactReference(
            "ePC-SAFT-validation@abc123:profiles/kij.json",
            "b" * 64,
        ),
        bootstrap_artifact=None,
        uncertainty_artifact=None,
        validation_campaign_artifact=None,
    )

    first = result.to_json_bytes(prepared=prepared, context=context)
    second = result.to_json_bytes(prepared=prepared, context=context)
    record = json.loads(first)

    assert first == second
    assert record["schema_id"] == "epcsaft-regression-result"
    assert record["schema_version"] == 1
    assert record["literature"]["reproduction_class"] == "MODERN_REFIT"
    assert record["literature"]["reconstruction_decisions"] == [
        "used the current exact-EOS scaled least-squares objective"
    ]
    assert record["problem"]["datasets"][0]["rows"][0]["acquisition"] == (
        "direct_measurement"
    )
    assert record["objective_identity"]["doi:example"]["row_weighting"] == (
        "observation_residual_scales"
    )
    assert record["parameters"][0]["final"] == result.parameters[0].final
    assert record["literature"]["source_printed_parameters"]["k_ij"] == "-0.01"
    assert record["installed_artifacts"]["distribution"] == "epcsaft==0.2.0.dev0"
    assert record["installed_artifacts"]["wheel_sha256"] == (
        "1567cda72e1b525526dc0e647af0c6fe711edcb70bc4cee08f06284e847956d9"
    )
    assert record["installed_artifacts"]["commit"] == (
        "7b97bab039e1c50a6f89522698af80493bea5f9e"
    )
    assert record["installed_artifacts"]["tree"] == (
        "d082a8f102b32705b6cd6669a3e31a8d4ea8acd0"
    )
    assert len(record["installed_artifacts"]["public_header_sha256"]) == 64
    assert len(record["installed_artifacts"]["static_library_sha256"]) == 64
    assert len(record["installed_artifacts"]["record_sha256"]) == 64
    assert record["status"]["authority_status"] == (
        "NO_AUTHORITY_CHANGE_REGRESSION_RESULT_ONLY"
    )
    assert record["literature"]["profile_artifact"]["artifact_id"] == (
        "ePC-SAFT-validation@abc123:profiles/kij.json"
    )
    assert record["literature"]["bootstrap_artifact"] is None
    assert record["literature"]["uncertainty_artifact"] is None
    assert record["literature"]["validation_campaign_artifact"] is None
    assert get_type_hints(RegressionResult)["problem"]
    assert get_type_hints(RegressionResult.to_record)["prepared"]
    assert get_type_hints(fit_parameters)["return"] is RegressionResult
    assert get_type_hints(type(prepared).evaluate)["return"] is PreparedFitEvaluation

    with pytest.raises(ValueError, match="finite"):
        replace(result, final_cost=float("nan")).to_json_bytes(prepared=prepared)
    with pytest.raises(ValueError, match="requires its exact PreparedFit"):
        result.to_record()
    with pytest.raises(ValueError, match="literature result export requires"):
        replace(result, preparation_fingerprint=None).to_record(context=context)
    tampered_dataset = ObservationDataset.from_records(
        FixedCompositionVleObservation,
        _records(),
        source=_source(),
        objective=_objective(),
        row_provenance={
            "row-1": replace(
                prepared.datasets[0].row_provenance[0][1],
                outlier_decision="tampered after fit",
            )
        },
    )
    with pytest.raises(ValueError, match="provenance does not own"):
        result.to_record(prepared=replace(prepared, datasets=(tampered_dataset,)))
    with pytest.raises(ValueError, match="complete model_identity"):
        ResultContext(reproduction_class=ReproductionClass.MODERN_REFIT)
    with pytest.raises(ValueError, match="reproduction class"):
        ResultContext(model_identity=context.model_identity)
    with pytest.raises(ValueError, match="source-printed parameters"):
        ResultContext(source_printed_parameters=(("k_ij", "-0.01"),))
    with pytest.raises(ValueError, match="explicit reconstruction decisions"):
        replace(context, reconstruction_decisions=())
    tuple_context = replace(
        context,
        reproduction_class=ReproductionClass.PUBLISHED_TUPLE_PROPERTY_REPLAY,
        reconstruction_decisions=(),
    )
    with pytest.raises(ValueError, match="no-fit property result"):
        result.to_record(prepared=prepared, context=tuple_context)
