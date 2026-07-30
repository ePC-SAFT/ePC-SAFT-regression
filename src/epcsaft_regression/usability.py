from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Iterable, Mapping
from dataclasses import asdict, dataclass, fields, is_dataclass
from enum import StrEnum
from types import MappingProxyType

from . import _native
from .parameter_regression import (
    AffineParameterTransform,
    AqueousKijMeanIonicActivityObservation,
    ComponentParameterIdentity,
    FixedCompositionVleObservation,
    IonSolvationKijObservation,
    MeanIonicActivityObservation,
    ModelParameterIdentity,
    ObservationPartition,
    PairParameterIdentity,
    ParameterCapability,
    ParameterCoordinate,
    ParameterFamily,
    PureDensityObservation,
    PureSaturationObservation,
    PureVaporPressureObservation,
    RegressionObservation,
    RegressionProblem,
    RegressionResult,
    RelativePermittivityRatioObservation,
    SolvationGibbsObservation,
    SourceDescriptor,
    UnsupportedParameterCapability,
    _evaluate_parameters,
    _require_finite,
    _require_nonempty_string,
    _require_sha256,
    canonical_dataset_sha256,
    fit_parameters,
    parameter_capabilities,
)


class AcquisitionClass(StrEnum):
    DIRECT_MEASUREMENT = "direct_measurement"
    DIGITIZED_VALUE = "digitized_value"
    DATABASE_RECORD = "database_record"
    AUTHOR_CORRELATION = "author_correlation"
    RECONSTRUCTED_CORRELATION = "reconstructed_correlation"


class ReproductionClass(StrEnum):
    EXACT_AUTHOR_METHOD_REPLAY = "EXACT_AUTHOR_METHOD_REPLAY"
    SOURCE_FAITHFUL_RECONSTRUCTION = "SOURCE_FAITHFUL_RECONSTRUCTION"
    PUBLISHED_TUPLE_PROPERTY_REPLAY = "PUBLISHED_TUPLE_PROPERTY_REPLAY"
    MODERN_REFIT = "MODERN_REFIT"


@dataclass(frozen=True, slots=True)
class CorrelationProvenance:
    equation: str
    coefficients: tuple[tuple[str, float], ...]
    units: str
    validity_interval: str
    sampling_grid: tuple[float, ...]
    transformation_record: str

    def __post_init__(self) -> None:
        for name in ("equation", "units", "validity_interval", "transformation_record"):
            _require_nonempty_string(getattr(self, name), name)
        if not self.coefficients or not self.sampling_grid:
            raise ValueError("correlation coefficients and sampling_grid must be nonempty")
        if any(
            type(value) not in (int, float) or not math.isfinite(value)
            for _, value in self.coefficients
        ) or any(
            type(value) not in (int, float) or not math.isfinite(value)
            for value in self.sampling_grid
        ):
            raise ValueError("correlation coefficients and sampling_grid must be finite")


@dataclass(frozen=True, slots=True)
class RowProvenance:
    acquisition: AcquisitionClass
    duplicate_decision: str
    exclusion_decision: str
    critical_region_decision: str
    censoring_decision: str
    outlier_decision: str
    correlation: CorrelationProvenance | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.acquisition, AcquisitionClass):
            raise TypeError("acquisition must be an AcquisitionClass")
        for name in (
            "duplicate_decision",
            "exclusion_decision",
            "critical_region_decision",
            "censoring_decision",
            "outlier_decision",
        ):
            _require_nonempty_string(getattr(self, name), name)
        correlation_backed = self.acquisition in (
            AcquisitionClass.AUTHOR_CORRELATION,
            AcquisitionClass.RECONSTRUCTED_CORRELATION,
        )
        if correlation_backed != (self.correlation is not None):
            raise ValueError(
                "correlation acquisition classes require correlation provenance, "
                "and other acquisition classes forbid it"
            )


@dataclass(frozen=True, slots=True)
class ObjectiveContract:
    residual_family: str
    interpretation: str
    row_weighting: str
    covariance_interpretation: str
    loss: str
    loss_parameters: tuple[tuple[str, float], ...]
    failed_row_treatment: str

    def __post_init__(self) -> None:
        for name in (
            "residual_family",
            "interpretation",
            "row_weighting",
            "covariance_interpretation",
            "failed_row_treatment",
        ):
            _require_nonempty_string(getattr(self, name), name)
        if self.loss != "squared":
            raise ValueError("the current Ceres problem supports only squared loss")
        if self.loss_parameters:
            raise ValueError("squared loss has no loss parameters")


@dataclass(frozen=True, slots=True)
class SourceInput:
    source_id: str
    citation: str
    durable_locator: str
    source_artifact_sha256: str
    transformation_record: str
    units_and_bases: str
    use_basis: str
    residual_scale_rationale: str

    def __post_init__(self) -> None:
        for name in (
            "source_id",
            "citation",
            "durable_locator",
            "transformation_record",
            "units_and_bases",
            "use_basis",
            "residual_scale_rationale",
        ):
            _require_nonempty_string(getattr(self, name), name)
        _require_sha256(self.source_artifact_sha256, "source_artifact_sha256")


def _json_value(value: object) -> object:
    if isinstance(value, StrEnum):
        return value.value
    if is_dataclass(value):
        return {key: _json_value(item) for key, item in asdict(value).items()}
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    return value


def _canonical_sha256(value: object) -> str:
    encoded = json.dumps(
        _json_value(value), sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode()
    return hashlib.sha256(encoded).hexdigest()


_TUPLE_FIELDS = {
    "component_ids",
    "active_pair_component_ids",
    "fixed_k_ij",
    "chemical_potential_scales",
    "liquid_volume_bounds_m3_per_mol",
    "vapor_volume_bounds_m3_per_mol",
    "volume_bounds_m3_per_mol",
}
_OBSERVATION_TYPES = (
    FixedCompositionVleObservation,
    PureSaturationObservation,
    PureVaporPressureObservation,
    PureDensityObservation,
    MeanIonicActivityObservation,
    AqueousKijMeanIonicActivityObservation,
    IonSolvationKijObservation,
    SolvationGibbsObservation,
    RelativePermittivityRatioObservation,
)


@dataclass(frozen=True, slots=True)
class ObservationDataset:
    source: SourceInput
    observations: tuple[RegressionObservation, ...]
    objective: ObjectiveContract
    row_provenance: tuple[tuple[str, RowProvenance], ...]
    provenance_sha256: str

    def __post_init__(self) -> None:
        if not self.observations:
            raise ValueError("dataset observations must be nonempty")
        row_ids = tuple(row.row_id for row in self.observations)
        if len(set(row_ids)) != len(row_ids):
            raise ValueError("dataset row IDs must be unique")
        provenance_ids = tuple(row_id for row_id, _ in self.row_provenance)
        if len(set(provenance_ids)) != len(provenance_ids) or set(
            provenance_ids
        ) != set(row_ids):
            raise ValueError("row provenance must match every observation row exactly")
        expected = _canonical_sha256(
            {
                "source": self.source,
                "objective": self.objective,
                "rows": tuple(
                    (row_id, provenance) for row_id, provenance in self.row_provenance
                ),
                "canonical_dataset_sha256": canonical_dataset_sha256(self.observations),
            }
        )
        if self.provenance_sha256 != expected:
            raise ValueError("provenance_sha256 does not match the validated dataset")

    @classmethod
    def from_records(
        cls,
        observation_type: type,
        records: Iterable[Mapping[str, object]],
        *,
        source: SourceInput,
        objective: ObjectiveContract,
        row_provenance: Mapping[str, RowProvenance],
    ) -> ObservationDataset:
        if observation_type not in _OBSERVATION_TYPES:
            raise ValueError("observation_type is not a fit-ready direct-EOS contract")
        allowed = {field.name for field in fields(observation_type)}
        observations: list[RegressionObservation] = []
        for index, record in enumerate(records):
            if not isinstance(record, Mapping):
                raise TypeError(f"row {index}: record must be a mapping")
            row_id = str(record.get("row_id", f"index-{index}"))
            extra = sorted(set(record) - allowed)
            missing = sorted(allowed - set(record))
            if extra or missing:
                detail = []
                if missing:
                    detail.append(f"missing {', '.join(missing)}")
                if extra:
                    detail.append(f"unexpected {', '.join(extra)}")
                raise ValueError(f"{row_id}: {'; '.join(detail)}")
            values = dict(record)
            for name in _TUPLE_FIELDS.intersection(values):
                if isinstance(values[name], list):
                    values[name] = tuple(values[name])
            if isinstance(values.get("partition"), str):
                try:
                    values["partition"] = ObservationPartition(values["partition"])
                except ValueError as exc:
                    raise ValueError(f"{row_id}: invalid partition") from exc
            try:
                observations.append(observation_type(**values))
            except (TypeError, ValueError) as exc:
                raise ValueError(f"{row_id}: {exc}") from exc
        rows = tuple(observations)
        mismatched_sources = tuple(
            row.row_id for row in rows if row.source_id != source.source_id
        )
        if mismatched_sources:
            raise ValueError(
                "rows reference a source_id other than the dataset source: "
                + ", ".join(mismatched_sources)
            )
        provenance = tuple(sorted(row_provenance.items()))
        payload = {
            "source": source,
            "objective": objective,
            "rows": provenance,
            "canonical_dataset_sha256": canonical_dataset_sha256(rows),
        }
        return cls(source, rows, objective, provenance, _canonical_sha256(payload))


@dataclass(frozen=True, slots=True)
class SolverControls:
    maximum_iterations: int
    maximum_solver_time_seconds: float
    function_tolerance: float
    gradient_tolerance: float
    parameter_tolerance: float

    def __post_init__(self) -> None:
        if type(self.maximum_iterations) is not int or self.maximum_iterations <= 0:
            raise ValueError("maximum_iterations must be positive")
        for name in (
            "maximum_solver_time_seconds",
            "function_tolerance",
            "gradient_tolerance",
            "parameter_tolerance",
        ):
            _require_finite(getattr(self, name), name, positive=True)


@dataclass(frozen=True, slots=True)
class RankControls:
    maximum_condition_number: float

    def __post_init__(self) -> None:
        _require_finite(
            self.maximum_condition_number,
            "maximum_condition_number",
            positive=True,
        )


@dataclass(frozen=True, slots=True)
class ConfirmationControls:
    parameter_scaled_max_delta: float
    cost_relative_delta: float

    def __post_init__(self) -> None:
        _require_finite(
            self.parameter_scaled_max_delta,
            "parameter_scaled_max_delta",
            positive=True,
        )
        _require_finite(
            self.cost_relative_delta,
            "cost_relative_delta",
            positive=True,
        )


@dataclass(frozen=True, slots=True)
class ParameterRequest:
    family: ParameterFamily
    identity: PairParameterIdentity | ComponentParameterIdentity | ModelParameterIdentity
    transform: AffineParameterTransform
    lower_bound: float
    upper_bound: float

    def __post_init__(self) -> None:
        if not isinstance(self.family, ParameterFamily):
            raise TypeError("family must be a ParameterFamily")
        if type(self.lower_bound) not in (int, float) or type(self.upper_bound) not in (
            int,
            float,
        ):
            raise ValueError("parameter bounds must be finite")
        if not (
            math.isfinite(self.lower_bound)
            and math.isfinite(self.upper_bound)
            and self.lower_bound < self.upper_bound
        ):
            raise ValueError("parameter bounds must be finite and increasing")


def support_view(
    model: object,
) -> tuple[ParameterCapability | UnsupportedParameterCapability, ...]:
    return parameter_capabilities(model)


def _unit(value: str) -> str:
    units = {"dimensionless": "1", "angstrom": "angstrom", "kelvin": "K"}
    if value not in units:
        raise ValueError(f"installed capability uses unsupported unit {value!r}")
    return units[value]


def _identity_matches(request: ParameterRequest, capability: ParameterCapability) -> bool:
    identity = request.identity.canonical_component_ids
    if isinstance(request.identity, ModelParameterIdentity):
        return capability.identity_shape == "model"
    if isinstance(request.identity, PairParameterIdentity):
        return identity == tuple(sorted(capability.active_component_ids))
    active = capability.active_component_ids or capability.component_ids
    return identity == active


_ASSOCIATING_FAMILIES = (
    ParameterFamily.SEGMENT_COUNT,
    ParameterFamily.SEGMENT_DIAMETER,
    ParameterFamily.DISPERSION_ENERGY_OVER_K,
    ParameterFamily.ASSOCIATION_ENERGY_OVER_K,
    ParameterFamily.ASSOCIATION_VOLUME,
)


def _resolve_coordinates(
    model: object, requests: tuple[ParameterRequest, ...]
) -> tuple[ParameterCoordinate, ...]:
    advertised = tuple(
        item
        for item in parameter_capabilities(model)
        if isinstance(item, ParameterCapability)
    )
    derivative_ready = tuple(
        item
        for item in advertised
        if item.installed_ready
    )
    if tuple(request.family for request in requests) == _ASSOCIATING_FAMILIES:
        matches = tuple(
            item
            for item in derivative_ready
            if item.capability_id == "neutral_pure_associating_joint_sigma_basis_v1"
        )
        if len(matches) != 1:
            raise ValueError("installed EOS does not advertise exactly one fixed-2B block")
        capability = matches[0]
        units = ("1", "angstrom", "K", "K", "1")
        return tuple(
            ParameterCoordinate(
                request.family,
                request.identity,
                capability.capability_id,
                capability.parameter_fingerprint,
                capability.topology_fingerprint,
                units[index],
                request.transform,
                request.lower_bound,
                request.upper_bound,
            )
            for index, request in enumerate(requests)
        )
    coordinates = []
    for request in requests:
        matches = tuple(
            item
            for item in derivative_ready
            if item.family is request.family and _identity_matches(request, item)
        )
        if len(matches) != 1:
            raise ValueError(
                f"installed EOS capability for {request.family.value} is "
                f"{'ambiguous' if len(matches) > 1 else 'unavailable'}"
            )
        capability = matches[0]
        coordinates.append(
            ParameterCoordinate(
                request.family,
                request.identity,
                capability.capability_id,
                capability.parameter_fingerprint,
                capability.topology_fingerprint,
                _unit(capability.coordinate_units[-1]),
                request.transform,
                request.lower_bound,
                request.upper_bound,
            )
        )
    return tuple(coordinates)


_OBJECTIVE_FAMILIES = {
    FixedCompositionVleObservation: "fixed_composition_vle",
    PureSaturationObservation: "pure_saturation",
    PureVaporPressureObservation: "pure_vapor_pressure",
    PureDensityObservation: "pure_density",
    MeanIonicActivityObservation: "mean_ionic_activity",
    AqueousKijMeanIonicActivityObservation: "aqueous_kij_mean_ionic_activity",
    IonSolvationKijObservation: "ion_solvation_kij",
    SolvationGibbsObservation: "solvation_gibbs",
    RelativePermittivityRatioObservation: "relative_permittivity_ratio",
}


def _validate_objective(dataset: ObservationDataset) -> None:
    families = {
        _OBJECTIVE_FAMILIES[type(row)] for row in dataset.observations
    }
    allowed = (
        families
        if len(families) == 1
        else {"pure_associating_mixed"}
        if families.issubset(
            {"pure_saturation", "pure_vapor_pressure", "pure_density"}
        )
        else set()
    )
    if dataset.objective.residual_family not in allowed:
        raise ValueError(
            "objective residual_family does not match the validated observation "
            f"contract; expected one of {sorted(allowed)}"
        )
    if dataset.objective.interpretation != "native_scaled_least_squares":
        raise ValueError(
            "objective interpretation must be 'native_scaled_least_squares'"
        )


@dataclass(frozen=True, slots=True)
class PreflightStart:
    start: tuple[float, ...]
    residual_count: int
    variable_count: int
    full_rank: int
    full_condition_number: float
    projected_parameter_rank: int
    projected_parameter_condition_number: float
    bounds_status: tuple[str, ...]
    derivative_complete: bool
    failure_reason: str


@dataclass(frozen=True, slots=True)
class PreflightReport:
    fitted_parameter_count: int
    lifted_variable_count: int
    residual_count: int
    training_row_count: int
    held_out_row_count: int
    stress_row_count: int
    partition_counts: tuple[tuple[str, int], ...]
    starts: tuple[PreflightStart, ...]
    local_evidence_only: bool
    practical_identifiability_gate: str | None
    ready: bool
    reasons: tuple[str, ...]


def _shape(rows: tuple[RegressionObservation, ...]) -> tuple[int, int]:
    residuals = lifted = 0
    for row in rows:
        if isinstance(
            row, (FixedCompositionVleObservation, PureSaturationObservation)
        ):
            residuals += 4
            lifted += 2
        elif isinstance(row, PureVaporPressureObservation):
            residuals += 3
            lifted += 2
        elif isinstance(row, PureDensityObservation):
            residuals += 2
            lifted += 1
        else:
            residuals += 1
    return residuals, lifted


def _lifted_start_variables(
    rows: tuple[RegressionObservation, ...],
) -> tuple[float, ...]:
    values: list[float] = []
    for row in rows:
        if isinstance(row, PureDensityObservation):
            values.append(
                math.log(
                    row.volume_start_m3_per_mol / row.volume_origin_m3_per_mol
                )
            )
        elif isinstance(
            row,
            (
                FixedCompositionVleObservation,
                PureSaturationObservation,
                PureVaporPressureObservation,
            ),
        ):
            values.extend(
                (
                    math.log(
                        row.liquid_volume_start_m3_per_mol
                        / row.liquid_volume_origin_m3_per_mol
                    ),
                    math.log(
                        row.vapor_volume_start_m3_per_mol
                        / row.vapor_volume_origin_m3_per_mol
                    ),
                )
            )
    return tuple(values)


@dataclass(frozen=True, slots=True)
class PreparedFit:
    model: object
    problem: RegressionProblem
    datasets: tuple[ObservationDataset, ...]

    def fit(self) -> RegressionResult:
        return fit_parameters(self.problem, self.model)

    def preflight(self) -> PreflightReport:
        parameter_count = len(self.problem.parameters)
        residual_count, lifted_count = _shape(self.problem.training_observations)
        lifted_start = _lifted_start_variables(self.problem.training_observations)
        if len(lifted_start) != lifted_count:
            raise ValueError("preflight lifted start accounting is incomplete")
        starts = []
        reasons: list[str] = []
        for physical_start, solver_start in zip(
            self.problem.start_vectors, self.problem.solver_start_vectors, strict=True
        ):
            bounds = tuple(
                "lower"
                if value == coordinate.lower_bound
                else "upper"
                if value == coordinate.upper_bound
                else "interior"
                for value, coordinate in zip(
                    physical_start, self.problem.parameters, strict=True
                )
            )
            try:
                residuals, jacobian = _evaluate_parameters(
                    self.problem,
                    self.model,
                    (*solver_start, *lifted_start),
                )
                variable_count = parameter_count + lifted_count
                complete = (
                    len(residuals) == residual_count
                    and len(jacobian) == residual_count * variable_count
                    and all(math.isfinite(value) for value in (*residuals, *jacobian))
                )
                if not complete:
                    raise ValueError("incomplete or nonfinite exact derivative payload")
                full, projected = _native.diagnose_jacobian(
                    parameter_count,
                    lifted_count,
                    residual_count,
                    jacobian,
                )
                full_rank, full_condition = int(full[1]), float(full[2])
                projected_rank, projected_condition = (
                    int(projected[1]),
                    float(projected[2]),
                )
                failure = ""
                if full_rank < variable_count:
                    failure = "rank_deficient_full_jacobian"
                elif projected_rank < parameter_count:
                    failure = "rank_deficient_projected_parameter_jacobian"
                elif (
                    full_condition > self.problem.maximum_condition_number
                    or projected_condition > self.problem.maximum_condition_number
                ):
                    failure = "poor_conditioning"
                if failure:
                    reasons.append(failure)
                starts.append(
                    PreflightStart(
                        physical_start,
                        residual_count,
                        variable_count,
                        full_rank,
                        full_condition,
                        projected_rank,
                        projected_condition,
                        bounds,
                        True,
                        failure,
                    )
                )
            except (ArithmeticError, RuntimeError, TypeError, ValueError) as exc:
                reason = f"eos_domain_or_derivative_failure: {exc}"
                reasons.append(reason)
                starts.append(
                    PreflightStart(
                        physical_start,
                        residual_count,
                        parameter_count + lifted_count,
                        0,
                        math.inf,
                        0,
                        math.inf,
                        bounds,
                        False,
                        reason,
                    )
                )
        gate = (
            "issue #28 nuisance-reoptimized profile and accepted-region evidence required"
            if tuple(parameter.family for parameter in self.problem.parameters)
            == _ASSOCIATING_FAMILIES
            else None
        )
        unique_reasons = tuple(dict.fromkeys(reasons))
        return PreflightReport(
            parameter_count,
            lifted_count,
            residual_count,
            len(self.problem.training_observations),
            len(self.problem.held_out_observations),
            len(self.problem.stress_observations),
            (
                ("training", len(self.problem.training_observations)),
                ("held_out", len(self.problem.held_out_observations)),
                ("stress", len(self.problem.stress_observations)),
            ),
            tuple(starts),
            True,
            gate,
            not unique_reasons,
            unique_reasons,
        )


def prepare_fit(
    model: object,
    *,
    datasets: tuple[ObservationDataset, ...],
    parameters: tuple[ParameterRequest, ...],
    parameter_slot_indices: tuple[int, ...],
    start_vectors: tuple[tuple[float, ...], ...],
    solver: SolverControls,
    rank: RankControls,
    confirmation: ConfirmationControls,
) -> PreparedFit:
    if not datasets:
        raise ValueError("datasets must be nonempty")
    for dataset in datasets:
        _validate_objective(dataset)
    coordinates = _resolve_coordinates(model, parameters)
    observations = tuple(row for dataset in datasets for row in dataset.observations)
    sources = tuple(
        SourceDescriptor(
            source_id=dataset.source.source_id,
            citation=dataset.source.citation,
            durable_locator=dataset.source.durable_locator,
            source_artifact_sha256=dataset.source.source_artifact_sha256,
            canonical_dataset_sha256=canonical_dataset_sha256(dataset.observations),
            transformation_record=dataset.source.transformation_record,
            units_and_bases=dataset.source.units_and_bases,
            use_basis=dataset.source.use_basis,
            residual_scale_rationale=dataset.source.residual_scale_rationale,
        )
        for dataset in datasets
    )
    problem = RegressionProblem(
        sources=sources,
        parameters=coordinates,
        parameter_slot_indices=parameter_slot_indices,
        start_vectors=start_vectors,
        observations=observations,
        maximum_condition_number=rank.maximum_condition_number,
        maximum_iterations=solver.maximum_iterations,
        maximum_solver_time_seconds=solver.maximum_solver_time_seconds,
        function_tolerance=solver.function_tolerance,
        gradient_tolerance=solver.gradient_tolerance,
        parameter_tolerance=solver.parameter_tolerance,
        confirmation_parameter_scaled_max_delta=confirmation.parameter_scaled_max_delta,
        confirmation_cost_relative_delta=confirmation.cost_relative_delta,
    )
    return PreparedFit(model, problem, datasets)


@dataclass(frozen=True, slots=True)
class ResultContext:
    reproduction_class: ReproductionClass | None = None
    model_identity: Mapping[str, object] | None = None
    source_printed_parameters: Mapping[str, str] | None = None
    practical_identifiability_artifact: Mapping[str, str] | None = None
    uncertainty_artifact: Mapping[str, str] | None = None
    data_identity: Mapping[str, object] | None = None
    objective_identity: Mapping[str, object] | None = None

    def __post_init__(self) -> None:
        if self.reproduction_class is not None and not isinstance(
            self.reproduction_class, ReproductionClass
        ):
            raise TypeError("reproduction_class must be a ReproductionClass")
        for name in (
            "model_identity",
            "source_printed_parameters",
            "practical_identifiability_artifact",
            "uncertainty_artifact",
            "data_identity",
            "objective_identity",
        ):
            value = getattr(self, name)
            if value is not None:
                object.__setattr__(self, name, MappingProxyType(dict(value)))


def _ensure_finite(value: object, path: str = "record") -> None:
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError(f"{path} must contain only finite numbers")
    if isinstance(value, Mapping):
        for key, item in value.items():
            _ensure_finite(item, f"{path}.{key}")
    elif isinstance(value, (tuple, list)):
        for index, item in enumerate(value):
            _ensure_finite(item, f"{path}[{index}]")


def _result_record(
    result: RegressionResult,
    *,
    prepared: PreparedFit | None = None,
    context: ResultContext | None = None,
) -> dict[str, object]:
    context = context or ResultContext()
    if prepared is not None and prepared.problem != result.problem:
        raise ValueError("prepared fit does not own this RegressionResult problem")
    positive = type(result.problem).__name__ == "PositiveEvaluatorProblem"
    if positive and (
        context.data_identity is None or context.objective_identity is None
    ):
        raise ValueError(
            "composed-positive export requires explicit data_identity and "
            "objective_identity"
        )
    problem = {
        "kind": type(result.problem).__name__,
        "resolved": _json_value(result.problem),
        "datasets": (
            [
                {
                    "source": _json_value(dataset.source),
                    "canonical_dataset_sha256": canonical_dataset_sha256(
                        dataset.observations
                    ),
                    "provenance_sha256": dataset.provenance_sha256,
                    "objective": _json_value(dataset.objective),
                    "rows": [
                        {
                            "row_id": row_id,
                            **_json_value(provenance),
                        }
                        for row_id, provenance in dataset.row_provenance
                    ],
                }
                for dataset in prepared.datasets
            ]
            if prepared is not None
            else None
        ),
    }
    record = {
        "schema_id": "epcsaft-regression-result",
        "schema_version": 1,
        "problem": problem,
        "data_identity": (
            _json_value(context.data_identity)
            if prepared is None
            else _json_value(
                {
                    dataset.source.source_id: {
                        row_id: provenance
                        for row_id, provenance in dataset.row_provenance
                    }
                    for dataset in prepared.datasets
                }
            )
        ),
        "objective_identity": (
            _json_value(context.objective_identity)
            if prepared is None
            else _json_value(
                {
                    dataset.source.source_id: dataset.objective
                    for dataset in prepared.datasets
                }
            )
        ),
        "capabilities": _json_value(result.capabilities),
        "provider_parameter_fingerprint": result.provider_parameter_fingerprint,
        "provider_topology_fingerprint": result.provider_topology_fingerprint,
        "status": {
            "solver_converged": result.solver_converged,
            "numerically_converged": result.numerically_converged,
            "workflow_valid": result.workflow_valid,
            "physical_status": result.physical_status,
            "scientific_status": result.scientific_status,
            "predictive_status": result.predictive_status,
            "termination": result.termination,
            "solution_usable": result.solution_usable,
            "failure_reasons": list(result.failure_reasons),
        },
        "solver": {
            "initial_cost": result.initial_cost,
            "final_cost": result.final_cost,
            "iterations": result.iterations,
            "residual_evaluation_count": result.residual_evaluation_count,
            "jacobian_evaluation_count": result.jacobian_evaluation_count,
        },
        "parameters": _json_value(result.parameters),
        "jacobian": _json_value(result.jacobian),
        "rows": _json_value(result.rows),
        "confirmation": {
            "count": result.confirmation_count,
            "parameter_scaled_max_delta": (
                result.confirmation_parameter_scaled_max_delta
            ),
            "cost_relative_max_delta": result.confirmation_cost_relative_max_delta,
            "usable": result.confirmations_usable,
        },
        "row_accounting": {
            "training": result.training_row_count,
            "held_out": result.held_out_row_count,
            "stress": result.stress_row_count,
            "evaluated": result.evaluated_row_count,
            "skipped": result.skipped_row_count,
            "failed": result.failed_row_count,
        },
        "literature": {
            "reproduction_class": (
                context.reproduction_class.value
                if context.reproduction_class is not None
                else None
            ),
            "model_identity": _json_value(context.model_identity),
            "source_printed_parameters": _json_value(
                context.source_printed_parameters
            ),
            "practical_identifiability_artifact": _json_value(
                context.practical_identifiability_artifact
            ),
            "uncertainty_artifact": _json_value(context.uncertainty_artifact),
        },
    }
    _ensure_finite(record)
    return record


def _result_json_bytes(
    result: RegressionResult,
    *,
    prepared: PreparedFit | None = None,
    context: ResultContext | None = None,
) -> bytes:
    return json.dumps(
        _result_record(result, prepared=prepared, context=context),
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")


__all__ = (
    "AcquisitionClass",
    "ConfirmationControls",
    "CorrelationProvenance",
    "ObjectiveContract",
    "ObservationDataset",
    "ParameterRequest",
    "PreflightReport",
    "PreflightStart",
    "PreparedFit",
    "RankControls",
    "ReproductionClass",
    "ResultContext",
    "RowProvenance",
    "SolverControls",
    "SourceInput",
    "prepare_fit",
    "support_view",
)
