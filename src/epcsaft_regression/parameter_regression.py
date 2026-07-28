from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
import hashlib
import json
import math
from typing import Iterable

from epcsaft import native_sdk

from . import _native


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


@dataclass(frozen=True, slots=True)
class ParameterCapability:
    capability_id: str
    family: ParameterFamily
    component_ids: tuple[str, ...]
    coordinate_kinds: tuple[str, ...]
    coordinate_units: tuple[str, ...]
    parameter_fingerprint: str
    topology_fingerprint: str
    derivative_order: int
    maturity: str
    authority_effect: str
    temperature_min_k: float
    temperature_max_k: float
    identity_shape: str
    observation_contract: str
    model_domain: str
    tensor_layout: str
    state_coordinate_count: int
    active_parameter_count: int
    helmholtz_basis_id: str
    unsupported_status: str
    domain_status: str
    active_component_ids: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class UnsupportedParameterCapability:
    capability_code: int
    schema_version: int
    parameter_family_code: int


@dataclass(frozen=True, slots=True)
class FittedParameterDiagnostic:
    family: ParameterFamily
    component_ids: tuple[str, ...]
    unit: str
    transform_origin: float
    transform_scale: float
    start: float
    final: float
    movement: float
    lower_bound: float
    upper_bound: float
    active_bound_distance: float
    active_bound: str | None


@dataclass(frozen=True, slots=True)
class GeneralJacobianDiagnostics:
    residual_count: int
    variable_count: int
    full_singular_values: tuple[float, ...]
    full_rank: int
    full_condition_number: float
    projected_parameter_singular_values: tuple[float, ...]
    projected_parameter_rank: int
    projected_parameter_condition_number: float


@dataclass(frozen=True, slots=True)
class GeneralRowDiagnostic:
    row_id: str
    partition: str
    liquid_volume_m3_per_mol: float
    vapor_volume_m3_per_mol: float
    scaled_residuals: tuple[float, float, float, float]
    observed_pressure_pa: float
    liquid_model_pressure_pa: float
    vapor_model_pressure_pa: float
    chemical_potential_differences_over_rt: tuple[float, float]
    derivative_status: str
    status: str
    evaluated: bool
    failure_reason: str


@dataclass(frozen=True, slots=True)
class PureSaturationRowDiagnostic:
    row_id: str
    partition: str
    liquid_volume_m3_per_mol: float
    vapor_volume_m3_per_mol: float
    scaled_residuals: tuple[float, float, float, float]
    observed_pressure_pa: float
    liquid_model_pressure_pa: float
    vapor_model_pressure_pa: float
    chemical_potential_difference_over_rt: float
    observed_liquid_density_kg_per_m3: float
    model_liquid_density_kg_per_m3: float
    derivative_status: str
    status: str
    evaluated: bool
    failure_reason: str


@dataclass(frozen=True, slots=True)
class DirectObservationRowDiagnostic:
    row_id: str
    partition: str
    observable: str
    observable_unit: str
    observed_value: float
    modeled_value: float
    scaled_residual: float
    provider_derivative: float
    derivative_status: str
    status: str
    evaluated: bool
    failure_reason: str


@dataclass(frozen=True, slots=True)
class RegressionResult:
    problem: RegressionProblem
    capability: ParameterCapability
    provider_parameter_fingerprint: str
    provider_topology_fingerprint: str
    solver_converged: bool
    numerically_converged: bool
    workflow_valid: bool
    physical_status: str
    scientific_status: str
    predictive_status: str
    termination: str
    solution_usable: bool
    initial_cost: float
    final_cost: float
    iterations: int
    parameter: FittedParameterDiagnostic
    jacobian: GeneralJacobianDiagnostics
    rows: tuple[
        GeneralRowDiagnostic
        | PureSaturationRowDiagnostic
        | DirectObservationRowDiagnostic,
        ...,
    ]
    confirmation_count: int
    confirmation_parameter_scaled_max_delta: float
    confirmation_cost_relative_max_delta: float
    confirmations_usable: bool
    training_row_count: int
    held_out_row_count: int
    stress_row_count: int
    evaluated_row_count: int
    skipped_row_count: int
    failed_row_count: int
    failure_reasons: tuple[str, ...]


def parameter_capabilities(
    model: object,
) -> tuple[ParameterCapability | UnsupportedParameterCapability, ...]:
    raw_capabilities = _native.parameter_capabilities(native_sdk(model))
    capabilities: list[ParameterCapability | UnsupportedParameterCapability] = []
    for raw in raw_capabilities:
        if len(raw) == 3 and type(raw[0]) is int:
            capabilities.append(
                UnsupportedParameterCapability(
                    capability_code=raw[0],
                    schema_version=raw[1],
                    parameter_family_code=raw[2],
                )
            )
            continue
        capabilities.append(
            ParameterCapability(
                capability_id=raw[0],
                family=ParameterFamily(raw[1]),
                component_ids=tuple(raw[2]),
                coordinate_kinds=tuple(raw[3]),
                coordinate_units=tuple(raw[4]),
                parameter_fingerprint=raw[5],
                topology_fingerprint=raw[6],
                derivative_order=raw[7],
                maturity=raw[8],
                authority_effect=raw[9],
                temperature_min_k=raw[10],
                temperature_max_k=raw[11],
                identity_shape=raw[12],
                observation_contract=raw[13],
                model_domain=raw[14],
                tensor_layout=raw[15],
                state_coordinate_count=raw[16],
                active_parameter_count=raw[17],
                helmholtz_basis_id=raw[18],
                unsupported_status=raw[19],
                domain_status=raw[20],
                active_component_ids=tuple(raw[21]),
            )
        )
    return tuple(capabilities)


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
class ComponentParameterIdentity:
    component_id: str

    def __post_init__(self) -> None:
        _require_nonempty_string(self.component_id, "component_id")

    @property
    def component_ids(self) -> tuple[str]:
        return (self.component_id,)

    @property
    def canonical_component_ids(self) -> tuple[str]:
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
    identity: PairParameterIdentity | ComponentParameterIdentity
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
        pair_families = (ParameterFamily.K_IJ, ParameterFamily.L_IJ)
        component_units = {
            ParameterFamily.SEGMENT_COUNT: "1",
            ParameterFamily.SEGMENT_DIAMETER: "angstrom",
            ParameterFamily.DISPERSION_ENERGY_OVER_K: "K",
            ParameterFamily.BORN_DIAMETER: "angstrom",
            ParameterFamily.SOLVATION_FACTOR: "1",
        }
        if self.family in pair_families:
            if not isinstance(self.identity, PairParameterIdentity):
                raise TypeError(
                    "pair-parameter identity must be a PairParameterIdentity"
                )
            expected_unit = "1"
        elif self.family in component_units:
            if not isinstance(self.identity, ComponentParameterIdentity):
                raise TypeError(
                    "component-parameter identity must be a "
                    "ComponentParameterIdentity"
                )
            expected_unit = component_units[self.family]
        else:
            raise ValueError(
                "the v1 coordinate contract supports only k_ij, l_ij, "
                "segment_count, segment_diameter, and "
                "dispersion_energy_over_k, born_diameter, or "
                "solvation_factor"
            )
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
        if self.unit != expected_unit:
            raise ValueError(
                f"{self.family.value} parameter unit must be {expected_unit!r}"
            )
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


@dataclass(frozen=True, slots=True)
class PureSaturationObservation:
    row_id: str
    source_id: str
    source_locator: str
    component_id: str
    temperature_k: float
    pressure_pa: float
    liquid_density_kg_per_m3: float
    molar_mass_kg_per_mol: float
    pressure_scale_pa: float
    chemical_potential_scale: float
    liquid_density_scale_kg_per_m3: float
    liquid_volume_origin_m3_per_mol: float
    liquid_volume_start_m3_per_mol: float
    liquid_volume_bounds_m3_per_mol: tuple[float, float]
    vapor_volume_origin_m3_per_mol: float
    vapor_volume_start_m3_per_mol: float
    vapor_volume_bounds_m3_per_mol: tuple[float, float]
    partition: ObservationPartition

    def __post_init__(self) -> None:
        for field in ("row_id", "source_id", "source_locator", "component_id"):
            _require_nonempty_string(getattr(self, field), field)
        for field in (
            "temperature_k",
            "pressure_pa",
            "liquid_density_kg_per_m3",
            "molar_mass_kg_per_mol",
            "pressure_scale_pa",
            "chemical_potential_scale",
            "liquid_density_scale_kg_per_m3",
        ):
            _require_finite(getattr(self, field), field, positive=True)
        self._validate_volume("liquid")
        self._validate_volume("vapor")
        if not isinstance(self.partition, ObservationPartition):
            raise TypeError("partition must be an ObservationPartition")

    @property
    def component_ids(self) -> tuple[str]:
        return (self.component_id,)

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


def _validate_direct_component_identity(
    component_ids: tuple[str, ...], active_component_id: str
) -> None:
    if type(component_ids) is not tuple or len(component_ids) < 2:
        raise ValueError("component_ids must contain the ordered Provider model")
    for component_id in component_ids:
        _require_nonempty_string(component_id, "component_id")
    if len(set(component_ids)) != len(component_ids):
        raise ValueError("component_ids must be distinct")
    _require_nonempty_string(active_component_id, "active_component_id")
    if active_component_id not in component_ids:
        raise ValueError("active_component_id must occur in component_ids")


@dataclass(frozen=True, slots=True)
class MeanIonicActivityObservation:
    row_id: str
    source_id: str
    source_locator: str
    component_ids: tuple[str, ...]
    active_component_id: str
    temperature_k: float
    pressure_pa: float
    formula_unit_molality_mol_per_kg: float
    observed_mean_ionic_activity_coefficient: float
    relative_residual_scale: float
    partition: ObservationPartition

    def __post_init__(self) -> None:
        for field in ("row_id", "source_id", "source_locator"):
            _require_nonempty_string(getattr(self, field), field)
        _validate_direct_component_identity(
            self.component_ids, self.active_component_id
        )
        for field in (
            "temperature_k",
            "pressure_pa",
            "formula_unit_molality_mol_per_kg",
            "observed_mean_ionic_activity_coefficient",
            "relative_residual_scale",
        ):
            _require_finite(getattr(self, field), field, positive=True)
        if not isinstance(self.partition, ObservationPartition):
            raise TypeError("partition must be an ObservationPartition")


@dataclass(frozen=True, slots=True)
class SolvationGibbsObservation:
    row_id: str
    source_id: str
    source_locator: str
    component_ids: tuple[str, ...]
    active_component_id: str
    temperature_k: float
    pressure_pa: float
    observed_solvation_gibbs_j_per_mol: float
    residual_scale_j_per_mol: float
    partition: ObservationPartition

    def __post_init__(self) -> None:
        for field in ("row_id", "source_id", "source_locator"):
            _require_nonempty_string(getattr(self, field), field)
        _validate_direct_component_identity(
            self.component_ids, self.active_component_id
        )
        _require_finite(self.temperature_k, "temperature_k", positive=True)
        _require_finite(self.pressure_pa, "pressure_pa", positive=True)
        _require_finite(
            self.observed_solvation_gibbs_j_per_mol,
            "observed_solvation_gibbs_j_per_mol",
        )
        _require_finite(
            self.residual_scale_j_per_mol,
            "residual_scale_j_per_mol",
            positive=True,
        )
        if not isinstance(self.partition, ObservationPartition):
            raise TypeError("partition must be an ObservationPartition")


DirectObservation = MeanIonicActivityObservation | SolvationGibbsObservation
RegressionObservation = (
    FixedCompositionVleObservation
    | PureSaturationObservation
    | DirectObservation
)


def _canonical_row(row: RegressionObservation) -> dict[str, object]:
    if isinstance(row, MeanIonicActivityObservation):
        return {
            "row_id": row.row_id,
            "source_id": row.source_id,
            "source_locator": row.source_locator,
            "component_ids": list(row.component_ids),
            "active_component_id": row.active_component_id,
            "temperature_k": row.temperature_k,
            "pressure_pa": row.pressure_pa,
            "formula_unit_molality_mol_per_kg": (
                row.formula_unit_molality_mol_per_kg
            ),
            "observed_mean_ionic_activity_coefficient": (
                row.observed_mean_ionic_activity_coefficient
            ),
            "relative_residual_scale": row.relative_residual_scale,
            "partition": row.partition.value,
        }
    if isinstance(row, SolvationGibbsObservation):
        return {
            "row_id": row.row_id,
            "source_id": row.source_id,
            "source_locator": row.source_locator,
            "component_ids": list(row.component_ids),
            "active_component_id": row.active_component_id,
            "temperature_k": row.temperature_k,
            "pressure_pa": row.pressure_pa,
            "observed_solvation_gibbs_j_per_mol": (
                row.observed_solvation_gibbs_j_per_mol
            ),
            "residual_scale_j_per_mol": row.residual_scale_j_per_mol,
            "partition": row.partition.value,
        }
    if isinstance(row, PureSaturationObservation):
        return {
            "row_id": row.row_id,
            "source_id": row.source_id,
            "source_locator": row.source_locator,
            "component_ids": list(row.component_ids),
            "temperature_k": row.temperature_k,
            "pressure_pa": row.pressure_pa,
            "liquid_density_kg_per_m3": row.liquid_density_kg_per_m3,
            "molar_mass_kg_per_mol": row.molar_mass_kg_per_mol,
            "pressure_scale_pa": row.pressure_scale_pa,
            "chemical_potential_scale": row.chemical_potential_scale,
            "liquid_density_scale_kg_per_m3": (
                row.liquid_density_scale_kg_per_m3
            ),
            "liquid_volume_origin_m3_per_mol": (
                row.liquid_volume_origin_m3_per_mol
            ),
            "liquid_volume_start_m3_per_mol": (
                row.liquid_volume_start_m3_per_mol
            ),
            "liquid_volume_bounds_m3_per_mol": list(
                row.liquid_volume_bounds_m3_per_mol
            ),
            "vapor_volume_origin_m3_per_mol": (
                row.vapor_volume_origin_m3_per_mol
            ),
            "vapor_volume_start_m3_per_mol": (
                row.vapor_volume_start_m3_per_mol
            ),
            "vapor_volume_bounds_m3_per_mol": list(
                row.vapor_volume_bounds_m3_per_mol
            ),
            "partition": row.partition.value,
        }
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
    observations: Iterable[RegressionObservation],
) -> str:
    rows = sorted((_canonical_row(row) for row in observations), key=lambda row: row["row_id"])
    payload = json.dumps(rows, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


@dataclass(frozen=True, slots=True)
class RegressionProblem:
    sources: tuple[SourceDescriptor, ...]
    parameters: tuple[ParameterCoordinate, ...]
    observations: tuple[RegressionObservation, ...]
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
            row_identity = (
                (row.active_component_id,)
                if isinstance(
                    row,
                    (MeanIonicActivityObservation, SolvationGibbsObservation),
                )
                else (
                    tuple(sorted(row.component_ids))
                    if len(row.component_ids) == 2
                    else row.component_ids
                )
            )
            if row_identity != identity:
                raise ValueError(
                    f"observation {row.row_id!r} identity does not match the parameter identity"
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
    def training_observations(self) -> tuple[RegressionObservation, ...]:
        return tuple(
            row for row in self.observations if row.partition is ObservationPartition.TRAINING
        )

    @property
    def held_out_observations(self) -> tuple[RegressionObservation, ...]:
        return tuple(
            row for row in self.observations if row.partition is ObservationPartition.HELD_OUT
        )

    @property
    def stress_observations(self) -> tuple[RegressionObservation, ...]:
        return tuple(
            row for row in self.observations if row.partition is ObservationPartition.STRESS
        )


def _row_payload(row: RegressionObservation) -> tuple[object, ...]:
    if isinstance(row, MeanIonicActivityObservation):
        return (
            row.row_id,
            row.partition.value,
            row.temperature_k,
            row.pressure_pa,
            row.formula_unit_molality_mol_per_kg,
            row.observed_mean_ionic_activity_coefficient,
            row.relative_residual_scale,
        )
    if isinstance(row, SolvationGibbsObservation):
        return (
            row.row_id,
            row.partition.value,
            row.temperature_k,
            row.pressure_pa,
            0.0,
            row.observed_solvation_gibbs_j_per_mol,
            row.residual_scale_j_per_mol,
        )
    if isinstance(row, PureSaturationObservation):
        return (
            row.row_id,
            row.partition.value,
            row.temperature_k,
            row.pressure_pa,
            row.pressure_scale_pa,
            row.chemical_potential_scale,
            row.molar_mass_kg_per_mol,
            row.liquid_density_kg_per_m3,
            row.liquid_density_scale_kg_per_m3,
            row.liquid_volume_origin_m3_per_mol,
            row.liquid_volume_start_m3_per_mol,
            row.liquid_volume_bounds_m3_per_mol[0],
            row.liquid_volume_bounds_m3_per_mol[1],
            row.vapor_volume_origin_m3_per_mol,
            row.vapor_volume_start_m3_per_mol,
            row.vapor_volume_bounds_m3_per_mol[0],
            row.vapor_volume_bounds_m3_per_mol[1],
        )
    return (
        row.row_id,
        row.partition.value,
        row.temperature_k,
        row.pressure_pa,
        row.liquid_mole_fraction_first,
        row.vapor_mole_fraction_first,
        row.pressure_scale_pa,
        row.chemical_potential_scales[0],
        row.chemical_potential_scales[1],
        row.liquid_volume_origin_m3_per_mol,
        row.liquid_volume_start_m3_per_mol,
        row.liquid_volume_bounds_m3_per_mol[0],
        row.liquid_volume_bounds_m3_per_mol[1],
        row.vapor_volume_origin_m3_per_mol,
        row.vapor_volume_start_m3_per_mol,
        row.vapor_volume_bounds_m3_per_mol[0],
        row.vapor_volume_bounds_m3_per_mol[1],
    )


def _native_payload(
    problem: RegressionProblem, capability: ParameterCapability
) -> tuple[object, ...]:
    parameter = problem.parameters[0]
    return (
        parameter.capability_id,
        parameter.provider_parameter_fingerprint,
        parameter.provider_topology_fingerprint,
        capability.component_ids,
        parameter.transform.origin,
        parameter.transform.scale,
        parameter.lower_bound,
        parameter.upper_bound,
        parameter.starts,
        problem.maximum_condition_number,
        problem.maximum_iterations,
        problem.function_tolerance,
        problem.gradient_tolerance,
        problem.parameter_tolerance,
        problem.confirmation_parameter_scaled_max_delta,
        problem.confirmation_cost_relative_delta,
        tuple(_row_payload(row) for row in problem.training_observations),
        tuple(
            _row_payload(row)
            for row in (*problem.held_out_observations, *problem.stress_observations)
        ),
    )


def _matched_capability(
    problem: RegressionProblem, model: object
) -> ParameterCapability:
    capabilities = parameter_capabilities(model)
    parameter = problem.parameters[0]
    matches = tuple(
        capability
        for capability in capabilities
        if isinstance(capability, ParameterCapability)
        and capability.capability_id == parameter.capability_id
    )
    if len(matches) != 1:
        raise ValueError(
            "installed Provider does not advertise the requested parameter capability"
        )
    capability = matches[0]
    if (
        capability.family is not parameter.family
        or capability.parameter_fingerprint
        != parameter.provider_parameter_fingerprint
        or capability.topology_fingerprint
        != parameter.provider_topology_fingerprint
    ):
        raise ValueError(
            "regression parameter does not match the installed Provider capability"
        )
    provider_unit = {
        "dimensionless": "1",
        "angstrom": "angstrom",
        "kelvin": "K",
    }.get(capability.coordinate_units[-1])
    if provider_unit != parameter.unit:
        raise ValueError(
            "regression parameter unit does not match the installed Provider "
            "capability"
        )
    for row in problem.observations:
        if row.component_ids != capability.component_ids:
            raise ValueError(
                f"observation {row.row_id!r} component order does not match the Provider model"
            )
        if not (
            capability.temperature_min_k
            <= row.temperature_k
            <= capability.temperature_max_k
        ):
            raise ValueError(
                f"observation {row.row_id!r} temperature is outside the Provider capability domain"
            )
        if isinstance(row, MeanIonicActivityObservation):
            if (
                capability.observation_contract
                != "aqueous_mean_ionic_activity"
                or capability.active_component_ids
                != parameter.identity.canonical_component_ids
            ):
                raise ValueError(
                    "mean-ionic-activity observation does not match the "
                    "Provider direct-observable capability"
                )
        elif isinstance(row, SolvationGibbsObservation):
            if (
                capability.observation_contract != "ion_solvation_gibbs"
                or capability.active_component_ids
                != parameter.identity.canonical_component_ids
            ):
                raise ValueError(
                    "solvation-Gibbs observation does not match the Provider "
                    "direct-observable capability"
                )
        elif capability.observation_contract != (
            "fixed_composition_helmholtz_phase"
        ):
            raise ValueError(
                "phase observation does not match the Provider capability"
            )
    return capability


def _evaluate_parameters(
    problem: RegressionProblem, model: object, variables: tuple[float, ...]
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    capability = _matched_capability(problem, model)
    residuals, jacobian = _native.evaluate_general(
        native_sdk(model), _native_payload(problem, capability), variables
    )
    return tuple(residuals), tuple(jacobian)


def fit_parameters(problem: RegressionProblem, model: object) -> RegressionResult:
    capability = _matched_capability(problem, model)
    native = _native.solve_general(
        native_sdk(model), _native_payload(problem, capability)
    )
    parameter_coordinate = problem.parameters[0]
    observations = {row.row_id: row for row in problem.observations}
    rows: tuple[
        GeneralRowDiagnostic
        | PureSaturationRowDiagnostic
        | DirectObservationRowDiagnostic,
        ...,
    ] = tuple(
        (
            DirectObservationRowDiagnostic(
                row_id=row[0],
                partition=row[1],
                observable=(
                    "mean_ionic_activity_coefficient"
                    if isinstance(
                        observation := observations[row[0]],
                        MeanIonicActivityObservation,
                    )
                    else "solvation_gibbs_energy"
                ),
                observable_unit=(
                    "1"
                    if isinstance(observation, MeanIonicActivityObservation)
                    else "J/mol"
                ),
                observed_value=(
                    observation.observed_mean_ionic_activity_coefficient
                    if isinstance(observation, MeanIonicActivityObservation)
                    else observation.observed_solvation_gibbs_j_per_mol
                ),
                modeled_value=row[7],
                scaled_residual=row[4][0],
                provider_derivative=row[8],
                derivative_status=(
                    "EXACT_PROVIDER_FIRST_DERIVATIVE"
                    if row[5]
                    else "UNAVAILABLE"
                ),
                status="evaluated" if row[5] else "failed",
                evaluated=row[5],
                failure_reason=row[6],
            )
            if isinstance(
                observations[row[0]],
                (MeanIonicActivityObservation, SolvationGibbsObservation),
            )
            else PureSaturationRowDiagnostic(
                row_id=row[0],
                partition=row[1],
                liquid_volume_m3_per_mol=row[2],
                vapor_volume_m3_per_mol=row[3],
                scaled_residuals=tuple(row[4]),
                observed_pressure_pa=observation.pressure_pa,
                liquid_model_pressure_pa=(
                    observation.pressure_pa
                    + row[4][0] * observation.pressure_scale_pa
                ),
                vapor_model_pressure_pa=(
                    observation.pressure_pa
                    + row[4][1] * observation.pressure_scale_pa
                ),
                chemical_potential_difference_over_rt=(
                    row[4][2] * observation.chemical_potential_scale
                ),
                observed_liquid_density_kg_per_m3=(
                    observation.liquid_density_kg_per_m3
                ),
                model_liquid_density_kg_per_m3=(
                    observation.molar_mass_kg_per_mol / row[2]
                ),
                derivative_status=(
                    "EXACT_PROVIDER_HESSIAN" if row[5] else "UNAVAILABLE"
                ),
                status="evaluated" if row[5] else "failed",
                evaluated=row[5],
                failure_reason=row[6],
            )
            if isinstance(
                observation := observations[row[0]],
                PureSaturationObservation,
            )
            else GeneralRowDiagnostic(
                row_id=row[0],
                partition=row[1],
                liquid_volume_m3_per_mol=row[2],
                vapor_volume_m3_per_mol=row[3],
                scaled_residuals=tuple(row[4]),
                observed_pressure_pa=observation.pressure_pa,
                liquid_model_pressure_pa=(
                    observation.pressure_pa
                    + row[4][0] * observation.pressure_scale_pa
                ),
                vapor_model_pressure_pa=(
                    observation.pressure_pa
                    + row[4][1] * observation.pressure_scale_pa
                ),
                chemical_potential_differences_over_rt=(
                    row[4][2] * observation.chemical_potential_scales[0],
                    row[4][3] * observation.chemical_potential_scales[1],
                ),
                derivative_status=(
                    "EXACT_PROVIDER_HESSIAN" if row[5] else "UNAVAILABLE"
                ),
                status="evaluated" if row[5] else "failed",
                evaluated=row[5],
                failure_reason=row[6],
            )
        )
        for row in native[20]
    )
    active_bound = native[7] or None
    parameter = FittedParameterDiagnostic(
        family=parameter_coordinate.family,
        component_ids=parameter_coordinate.identity.canonical_component_ids,
        unit=parameter_coordinate.unit,
        transform_origin=parameter_coordinate.transform.origin,
        transform_scale=parameter_coordinate.transform.scale,
        start=parameter_coordinate.starts[0],
        final=native[5],
        movement=native[5] - parameter_coordinate.starts[0],
        lower_bound=parameter_coordinate.lower_bound,
        upper_bound=parameter_coordinate.upper_bound,
        active_bound_distance=native[6],
        active_bound=active_bound,
    )
    jacobian = GeneralJacobianDiagnostics(
        residual_count=native[23],
        variable_count=native[22],
        full_singular_values=tuple(native[10]),
        full_rank=native[11],
        full_condition_number=native[12],
        projected_parameter_singular_values=tuple(native[13]),
        projected_parameter_rank=native[14],
        projected_parameter_condition_number=native[15],
    )
    solver_converged = (
        native[0] == "CONVERGENCE" and bool(native[1]) and not native[21]
    )
    confirmations_usable = bool(native[19])
    numerically_converged = (
        solver_converged
        and confirmations_usable
        and native[17] <= problem.confirmation_parameter_scaled_max_delta
        and native[18] <= problem.confirmation_cost_relative_delta
        and jacobian.projected_parameter_rank == len(problem.parameters)
        and math.isfinite(jacobian.full_condition_number)
        and jacobian.full_condition_number <= problem.maximum_condition_number
        and math.isfinite(jacobian.projected_parameter_condition_number)
        and jacobian.projected_parameter_condition_number
        <= problem.maximum_condition_number
    )
    workflow_valid = (
        all(row.evaluated and not row.failure_reason for row in rows)
        and len(rows) == len(problem.observations)
    )
    failures = tuple(
        reason
        for reason in (
            native[21],
            *(
                f"{row.row_id}: {row.failure_reason}"
                for row in rows
                if row.failure_reason
            ),
        )
        if reason
    )
    return RegressionResult(
        problem=problem,
        capability=capability,
        provider_parameter_fingerprint=capability.parameter_fingerprint,
        provider_topology_fingerprint=capability.topology_fingerprint,
        solver_converged=solver_converged,
        numerically_converged=numerically_converged,
        workflow_valid=workflow_valid,
        physical_status="NOT_ADJUDICATED_NO_ROW_ACCEPTANCE_CRITERIA",
        scientific_status="NOT_ADJUDICATED_NO_APPROVED_SCIENTIFIC_CUTOFF",
        predictive_status="NOT_ADJUDICATED_NO_APPROVED_HELD_OUT_CUTOFF",
        termination=native[0],
        solution_usable=bool(native[1]),
        initial_cost=native[2],
        final_cost=native[3],
        iterations=native[4],
        parameter=parameter,
        jacobian=jacobian,
        rows=rows,
        confirmation_count=native[16],
        confirmation_parameter_scaled_max_delta=native[17],
        confirmation_cost_relative_max_delta=native[18],
        confirmations_usable=confirmations_usable,
        training_row_count=len(problem.training_observations),
        held_out_row_count=len(problem.held_out_observations),
        stress_row_count=len(problem.stress_observations),
        evaluated_row_count=sum(row.evaluated for row in rows),
        skipped_row_count=0,
        failed_row_count=sum(not row.evaluated for row in rows),
        failure_reasons=failures,
    )
