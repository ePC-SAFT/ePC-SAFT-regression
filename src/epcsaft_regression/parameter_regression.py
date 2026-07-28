from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
import hashlib
import json
import math
from typing import Iterable


class ParameterFamily(StrEnum):
    SEGMENT_COUNT = "segment_count"
    SEGMENT_DIAMETER = "segment_diameter"
    DISPERSION_ENERGY_OVER_K = "dispersion_energy_over_k"
    RELATIVE_PERMITTIVITY = "relative_permittivity"
    BORN_DIAMETER = "born_diameter"
    SOLVATION_FACTOR = "solvation_factor"
    SCHRECKENBERG_DIELECTRIC_VOLUME = "schreckenberg_dielectric_volume"
    SCHRECKENBERG_DIELECTRIC_TEMPERATURE = "schreckenberg_dielectric_temperature"
    ZUBER_ION_SUPPRESSION_COEFFICIENT = "zuber_ion_suppression_coefficient"
    RUEBEN_DIPOLE_SCALING = "rueben_dipole_scaling"
    RUEBEN_POLARIZABILITY_SCALING = "rueben_polarizability_scaling"
    RUEBEN_CORRELATION_INTEGRAL_PARAMETER = "rueben_correlation_integral_parameter"
    K_IJ = "k_ij"
    L_IJ = "l_ij"
    ASSOCIATION_ENERGY_OVER_K = "association_energy_over_k"
    ASSOCIATION_VOLUME = "association_volume"
    DIELECTRIC_ION_SUPPRESSION_COEFFICIENT = "dielectric_ion_suppression_coefficient"
    IONIC_REGION_RELATIVE_PERMITTIVITY = "ionic_region_relative_permittivity"


class ObservationPartition(StrEnum):
    TRAINING = "training"
    HELD_OUT = "held_out"
    STRESS = "stress"


def _require_nonempty_string(value: str, field: str) -> None:
    if type(value) is not str or not value.strip():
        raise ValueError(f"{field} must be a nonempty string")


def _require_finite(value: float, field: str, *, positive: bool = False) -> None:
    if type(value) not in (int, float) or not math.isfinite(value):
        raise ValueError(f"{field} must be finite")
    if positive and value <= 0.0:
        raise ValueError(f"{field} must be positive")


def _require_sha256(value: str, field: str, *, prefixed: bool = False) -> None:
    _require_nonempty_string(value, field)
    prefix = "sha256:" if prefixed else ""
    body = value.removeprefix(prefix)
    if prefixed and not value.startswith(prefix):
        raise ValueError(f"{field} must begin with 'sha256:'")
    if len(body) != 64 or any(character not in "0123456789abcdef" for character in body):
        raise ValueError(f"{field} must be a lowercase SHA-256 identity")


@dataclass(frozen=True, slots=True)
class PairParameterIdentity:
    component_id_a: str
    component_id_b: str

    def __post_init__(self) -> None:
        _require_nonempty_string(self.component_id_a, "component_id_a")
        _require_nonempty_string(self.component_id_b, "component_id_b")
        if self.component_id_a == self.component_id_b:
            raise ValueError("pair parameter components must be distinct")
        canonical = tuple(sorted((self.component_id_a, self.component_id_b)))
        object.__setattr__(self, "component_id_a", canonical[0])
        object.__setattr__(self, "component_id_b", canonical[1])

    @property
    def component_ids(self) -> tuple[str, str]:
        return (self.component_id_a, self.component_id_b)

    @property
    def canonical_component_ids(self) -> tuple[str, str]:
        return self.component_ids


@dataclass(frozen=True, slots=True)
class AffineParameterTransform:
    origin: float
    scale: float

    def __post_init__(self) -> None:
        _require_finite(self.origin, "transform origin")
        _require_finite(self.scale, "transform scale")
        if self.scale == 0.0:
            raise ValueError("transform scale must be nonzero")

    def to_solver(self, physical_value: float) -> float:
        _require_finite(physical_value, "physical parameter value")
        return (physical_value - self.origin) / self.scale

    def to_physical(self, solver_value: float) -> float:
        _require_finite(solver_value, "solver parameter value")
        return self.origin + self.scale * solver_value


@dataclass(frozen=True, slots=True)
class ParameterCoordinate:
    family: ParameterFamily
    identity: PairParameterIdentity
    capability_id: str
    provider_parameter_fingerprint: str
    provider_topology_fingerprint: str
    unit: str
    transform: AffineParameterTransform
    lower_bound: float
    upper_bound: float
    starts: tuple[float, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.family, ParameterFamily):
            raise TypeError("family must be a ParameterFamily")
        if self.family is not ParameterFamily.K_IJ:
            raise ValueError("the v1 pair-coordinate contract currently supports only k_ij")
        if not isinstance(self.identity, PairParameterIdentity):
            raise TypeError("k_ij identity must be a PairParameterIdentity")
        _require_nonempty_string(self.capability_id, "capability_id")
        _require_sha256(
            self.provider_parameter_fingerprint,
            "provider parameter fingerprint",
            prefixed=True,
        )
        _require_sha256(
            self.provider_topology_fingerprint,
            "provider topology fingerprint",
            prefixed=True,
        )
        _require_nonempty_string(self.unit, "unit")
        if self.unit != "1":
            raise ValueError("k_ij unit must be '1'")
        if not isinstance(self.transform, AffineParameterTransform):
            raise TypeError("transform must be an AffineParameterTransform")
        _require_finite(self.lower_bound, "lower bound")
        _require_finite(self.upper_bound, "upper bound")
        if self.lower_bound >= self.upper_bound:
            raise ValueError("parameter bounds must be strictly increasing")
        if type(self.starts) is not tuple or len(self.starts) < 2:
            raise ValueError("starts must contain a primary and at least one confirmation start")
        for start in self.starts:
            _require_finite(start, "parameter start")
            if not self.lower_bound <= start <= self.upper_bound:
                raise ValueError("every parameter start must lie within the declared bounds")

    @property
    def solver_starts(self) -> tuple[float, ...]:
        return tuple(self.transform.to_solver(start) for start in self.starts)


@dataclass(frozen=True, slots=True)
class SourceDescriptor:
    source_id: str
    citation: str
    durable_locator: str
    source_artifact_sha256: str
    canonical_dataset_sha256: str
    transformation_record: str
    units_and_bases: str
    use_basis: str
    residual_scale_rationale: str

    def __post_init__(self) -> None:
        for field in (
            "source_id",
            "citation",
            "durable_locator",
            "transformation_record",
            "units_and_bases",
            "use_basis",
            "residual_scale_rationale",
        ):
            _require_nonempty_string(getattr(self, field), field)
        _require_sha256(self.source_artifact_sha256, "source artifact SHA-256")
        _require_sha256(self.canonical_dataset_sha256, "canonical dataset SHA-256")


@dataclass(frozen=True, slots=True)
class FixedCompositionVleObservation:
    row_id: str
    source_id: str
    source_locator: str
    component_ids: tuple[str, str]
    temperature_k: float
    pressure_pa: float
    liquid_mole_fraction_first: float
    vapor_mole_fraction_first: float
    pressure_scale_pa: float
    chemical_potential_scales: tuple[float, float]
    liquid_volume_origin_m3_per_mol: float
    liquid_volume_start_m3_per_mol: float
    liquid_volume_bounds_m3_per_mol: tuple[float, float]
    vapor_volume_origin_m3_per_mol: float
    vapor_volume_start_m3_per_mol: float
    vapor_volume_bounds_m3_per_mol: tuple[float, float]
    partition: ObservationPartition

    def __post_init__(self) -> None:
        for field in ("row_id", "source_id", "source_locator"):
            _require_nonempty_string(getattr(self, field), field)
        if type(self.component_ids) is not tuple or len(self.component_ids) != 2:
            raise ValueError("component_ids must contain exactly two component identifiers")
        for component_id in self.component_ids:
            _require_nonempty_string(component_id, "component_id")
        if self.component_ids[0] == self.component_ids[1]:
            raise ValueError("component_ids must be distinct")
        _require_finite(self.temperature_k, "temperature_k", positive=True)
        _require_finite(self.pressure_pa, "pressure_pa", positive=True)
        for value, field in (
            (self.liquid_mole_fraction_first, "liquid mole fraction"),
            (self.vapor_mole_fraction_first, "vapor mole fraction"),
        ):
            _require_finite(value, field)
            if not 0.0 < value < 1.0:
                raise ValueError(f"{field} must be strictly between zero and one")
        _require_finite(self.pressure_scale_pa, "pressure_scale_pa", positive=True)
        if type(self.chemical_potential_scales) is not tuple or len(
            self.chemical_potential_scales
        ) != 2:
            raise ValueError("chemical_potential_scales must contain two values")
        for scale in self.chemical_potential_scales:
            _require_finite(scale, "chemical potential scale", positive=True)
        self._validate_volume("liquid")
        self._validate_volume("vapor")
        if not isinstance(self.partition, ObservationPartition):
            raise TypeError("partition must be an ObservationPartition")

    def _validate_volume(self, phase: str) -> None:
        origin = getattr(self, f"{phase}_volume_origin_m3_per_mol")
        start = getattr(self, f"{phase}_volume_start_m3_per_mol")
        bounds = getattr(self, f"{phase}_volume_bounds_m3_per_mol")
        _require_finite(origin, f"{phase} volume origin", positive=True)
        _require_finite(start, f"{phase} volume start", positive=True)
        if type(bounds) is not tuple or len(bounds) != 2:
            raise ValueError(f"{phase} volume bounds must contain two values")
        lower, upper = bounds
        _require_finite(lower, f"{phase} volume lower bound", positive=True)
        _require_finite(upper, f"{phase} volume upper bound", positive=True)
        if lower >= upper:
            raise ValueError(f"{phase} volume bounds must be strictly increasing")
        if not lower <= start <= upper:
            raise ValueError(f"{phase} volume start must lie within its bounds")


def _canonical_row(row: FixedCompositionVleObservation) -> dict[str, object]:
    return {
        "row_id": row.row_id,
        "source_id": row.source_id,
        "source_locator": row.source_locator,
        "component_ids": list(row.component_ids),
        "temperature_k": row.temperature_k,
        "pressure_pa": row.pressure_pa,
        "liquid_mole_fraction_first": row.liquid_mole_fraction_first,
        "vapor_mole_fraction_first": row.vapor_mole_fraction_first,
        "pressure_scale_pa": row.pressure_scale_pa,
        "chemical_potential_scales": list(row.chemical_potential_scales),
        "liquid_volume_origin_m3_per_mol": row.liquid_volume_origin_m3_per_mol,
        "liquid_volume_start_m3_per_mol": row.liquid_volume_start_m3_per_mol,
        "liquid_volume_bounds_m3_per_mol": list(row.liquid_volume_bounds_m3_per_mol),
        "vapor_volume_origin_m3_per_mol": row.vapor_volume_origin_m3_per_mol,
        "vapor_volume_start_m3_per_mol": row.vapor_volume_start_m3_per_mol,
        "vapor_volume_bounds_m3_per_mol": list(row.vapor_volume_bounds_m3_per_mol),
        "partition": row.partition.value,
    }


def canonical_dataset_sha256(
    observations: Iterable[FixedCompositionVleObservation],
) -> str:
    rows = sorted((_canonical_row(row) for row in observations), key=lambda row: row["row_id"])
    payload = json.dumps(rows, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


@dataclass(frozen=True, slots=True)
class RegressionProblem:
    sources: tuple[SourceDescriptor, ...]
    parameters: tuple[ParameterCoordinate, ...]
    observations: tuple[FixedCompositionVleObservation, ...]
    maximum_condition_number: float
    maximum_iterations: int
    function_tolerance: float
    gradient_tolerance: float
    parameter_tolerance: float
    confirmation_parameter_scaled_max_delta: float
    confirmation_cost_relative_delta: float

    def __post_init__(self) -> None:
        if not self.sources:
            raise ValueError("at least one source is required")
        if not self.parameters:
            raise ValueError("at least one parameter is required")
        if not self.observations:
            raise ValueError("at least one observation is required")
        source_ids = tuple(source.source_id for source in self.sources)
        if len(set(source_ids)) != len(source_ids):
            raise ValueError("source_id values must be unique")
        row_ids = tuple(row.row_id for row in self.observations)
        if len(set(row_ids)) != len(row_ids):
            raise ValueError("duplicate row_id values are forbidden")
        parameter_keys = tuple(
            (parameter.family, parameter.identity.canonical_component_ids)
            for parameter in self.parameters
        )
        if len(set(parameter_keys)) != len(parameter_keys):
            raise ValueError("duplicate parameter identity values are forbidden")
        source_map = {source.source_id: source for source in self.sources}
        for row in self.observations:
            if row.source_id not in source_map:
                raise ValueError(f"observation {row.row_id!r} references an unknown source_id")
        for source in self.sources:
            rows = tuple(row for row in self.observations if row.source_id == source.source_id)
            if not rows:
                raise ValueError(f"source {source.source_id!r} has no observations")
            if canonical_dataset_sha256(rows) != source.canonical_dataset_sha256:
                raise ValueError(
                    f"source {source.source_id!r} canonical dataset SHA-256 does not match its rows"
                )
        if not self.training_observations:
            raise ValueError("at least one training observation is required")
        if len(self.parameters) != 1:
            raise ValueError("the first implementation accepts exactly one shared parameter")
        identity = self.parameters[0].identity.canonical_component_ids
        for row in self.observations:
            if tuple(sorted(row.component_ids)) != identity:
                raise ValueError(
                    f"observation {row.row_id!r} component pair does not match the parameter identity"
                )
        _require_finite(
            self.maximum_condition_number,
            "maximum_condition_number",
            positive=True,
        )
        if type(self.maximum_iterations) is not int or self.maximum_iterations <= 0:
            raise ValueError("maximum_iterations must be a positive integer")
        for field in (
            "function_tolerance",
            "gradient_tolerance",
            "parameter_tolerance",
            "confirmation_parameter_scaled_max_delta",
            "confirmation_cost_relative_delta",
        ):
            _require_finite(getattr(self, field), field, positive=True)

    @property
    def training_observations(self) -> tuple[FixedCompositionVleObservation, ...]:
        return tuple(
            row for row in self.observations if row.partition is ObservationPartition.TRAINING
        )

    @property
    def held_out_observations(self) -> tuple[FixedCompositionVleObservation, ...]:
        return tuple(
            row for row in self.observations if row.partition is ObservationPartition.HELD_OUT
        )

    @property
    def stress_observations(self) -> tuple[FixedCompositionVleObservation, ...]:
        return tuple(
            row for row in self.observations if row.partition is ObservationPartition.STRESS
        )
