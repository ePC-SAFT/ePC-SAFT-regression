from __future__ import annotations

from dataclasses import dataclass, replace
import math

from epcsaft import EPCSAFT, ParameterBundle, native_sdk, unit_registry
from epcsaft.records import (
    AssociationParameterRecord,
    ConstantCorrelation,
    ConstantPlusSumOfExponentialsCorrelation,
    ModelParameterRecord,
    PairParameterRecord,
    SingleParameterRecord,
    SiteRecord,
)

from . import _native
from .records import (
    FIGIEL_BORN_DIAMETER_TRACER_V1,
    FIGIEL_AQUEOUS_COMPONENTS,
    FIGIEL_AQUEOUS_KIJ_COORDINATES,
    FIGIEL_AQUEOUS_PUBLISHED_KIJ,
    FIGIEL_AQUEOUS_SALTS,
    FIGIEL_STAGED_AQUEOUS_RECOVERY_V1,
    BornDiameterTracerSpecification,
    FigielStagedAqueousRecoverySpecification,
    PureSaturationDataset,
    PureSaturationFitSpecification,
    SaturationObservation,
)


PROVIDER_CAPSULE = "epcsaft.native_sdk.v1"
PARAMETER_TRANSFORM = "p_j = start_j + parameter_scale_j * z_j"
LIQUID_VOLUME_TRANSFORM = "V_liquid = (molar_mass / observed_liquid_density) * exp(u_liquid)"
VAPOR_VOLUME_TRANSFORM = "V_vapor = (R * T / observed_pressure) * exp(u_vapor)"
REPORTING_PRESSURE_TRANSFORM = "P_report = observed_pressure * exp(u_pressure)"
PREDICTIVE_STATUS = "NOT_ADJUDICATED_NO_APPROVED_HELD_OUT_CUTOFF"
DIAMETER_TRANSFORM = "d_i = 3.0 angstrom + 1.0 angstrom * z_i"
BORN_RESIDUAL = "r_i = (G_i(d_i) - G_i_target) / abs(G_i_target)"
BORN_JACOBIAN = "J_ij = delta_ij * G_i_prime(d_i) * 1 angstrom / abs(G_i_target)"
AQUEOUS_MIAC_RESIDUAL = "r_q = 1 - gamma_q_model / gamma_q_observed"
AQUEOUS_MIAC_JACOBIAN = (
    "dr_q/dtheta_j = -(gamma_q_model/gamma_q_observed) "
    "* dln(gamma_q_model)/dtheta_j"
)
AQUEOUS_KIJ_COLUMNS = {
    "LiCl": (0, 3, 5),
    "NaCl": (1, 3, 6),
    "KCl": (2, 3, 7),
    "LiBr": (0, 4, 8),
    "NaBr": (1, 4, 9),
    "KBr": (2, 4, 10),
}


def _row_payload(row: SaturationObservation) -> tuple[object, ...]:
    return (
        row.row_id,
        row.component_id,
        row.temperature_k,
        row.pressure_pa,
        row.liquid_density_kg_m3,
        row.source_id,
    )


@dataclass(frozen=True, slots=True)
class ParameterDiagnostic:
    name: str
    unit: str
    start: float
    final: float
    movement: float
    lower_bound: float
    upper_bound: float
    active_bound: str | None


@dataclass(frozen=True, slots=True)
class JacobianDiagnostics:
    complete_columns: bool
    full_singular_values: tuple[float, ...]
    full_rank: int
    full_condition_number: float
    parameter_singular_values: tuple[float, ...]
    parameter_rank: int
    parameter_condition_number: float


@dataclass(frozen=True, slots=True)
class TrainingRowDiagnostic:
    row_id: str
    temperature_k: float
    observed_pressure_pa: float
    observed_liquid_density_kg_m3: float
    liquid_volume_m3: float
    vapor_volume_m3: float
    liquid_molar_density_mol_m3: float
    vapor_molar_density_mol_m3: float
    liquid_mass_density_kg_m3: float
    vapor_mass_density_kg_m3: float
    liquid_pressure_pa: float
    vapor_pressure_pa: float
    liquid_chemical_potential_over_rt: float
    vapor_chemical_potential_over_rt: float
    liquid_stability_slope: float
    vapor_stability_slope: float
    raw_residuals: tuple[float, float, float, float]
    scaled_residuals: tuple[float, float, float, float]


@dataclass(frozen=True, slots=True)
class ReportingRowDiagnostic:
    row_id: str
    temperature_k: float
    training: bool
    partition: str
    observed_pressure_pa: float
    predicted_pressure_pa: float
    pressure_relative_error: float
    observed_liquid_density_kg_m3: float
    predicted_liquid_density_kg_m3: float
    liquid_density_relative_error: float
    liquid_volume_m3: float
    vapor_volume_m3: float
    liquid_molar_density_mol_m3: float
    vapor_molar_density_mol_m3: float
    liquid_mass_density_kg_m3: float
    vapor_mass_density_kg_m3: float
    liquid_stability_slope: float
    vapor_stability_slope: float
    raw_equilibrium_residuals: tuple[float, float, float]
    termination: str
    solution_usable: bool
    physically_valid: bool
    failure_reasons: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class PureSaturationFitResult:
    component_id: str
    dataset_id: str
    specification_id: str
    provider_fingerprint: str
    compiled_problem_identity: tuple[str, ...]
    solver_converged: bool
    numerically_converged: bool
    physically_valid: bool
    predictive_status: str
    termination: str
    solution_usable: bool
    initial_cost: float
    final_cost: float
    iterations: int
    parameters: tuple[ParameterDiagnostic, ParameterDiagnostic, ParameterDiagnostic]
    jacobian: JacobianDiagnostics
    training_rows: tuple[TrainingRowDiagnostic, ...]
    reporting_rows: tuple[ReportingRowDiagnostic, ...]
    confirmation_termination: str
    confirmation_solution_usable: bool
    confirmation_parameter_scaled_max_delta: float
    confirmation_cost_relative_delta: float
    failure_reasons: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class BornObservationDiagnostic:
    target_id: str
    ion_label: str
    target_j_per_mol: float
    modeled_j_per_mol: float
    raw_error_j_per_mol: float
    scaled_residual: float
    derivative_j_per_mol_per_angstrom: float
    scaled_jacobian: float
    reference_molality_mol_per_kg: float
    reference_convergence_error: float
    provider_fingerprint: str


@dataclass(frozen=True, slots=True)
class BornStartDiagnostic:
    name: str
    termination: str
    solution_usable: bool
    initial_cost: float
    final_cost: float
    iterations: int
    transformed_parameters: tuple[float, ...]
    final_diameters_angstrom: tuple[float, ...]
    observations: tuple[BornObservationDiagnostic, ...]
    singular_values: tuple[float, ...]
    rank_threshold: float
    rank: int
    condition_number: float
    complete_columns: bool
    inactive_bounds: bool
    solver_converged: bool
    failure_reasons: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class BornParameterDiagnostic:
    ion_label: str
    active_component_id: str
    final_diameter_angstrom: float
    published_diameter_angstrom: float
    published_delta_angstrom: float
    lower_bound_angstrom: float
    upper_bound_angstrom: float
    scaled_lower_bound_distance: float
    scaled_upper_bound_distance: float
    active_bound: bool


@dataclass(frozen=True, slots=True)
class BornDiameterFitResult:
    specification_id: str
    compiled_problem_identity: tuple[str, ...]
    provider_fingerprints: tuple[str, ...]
    solver_converged: bool
    numerically_converged: bool
    workflow_valid: bool
    scientifically_valid: bool
    predictive_status: str
    parameters: tuple[BornParameterDiagnostic, ...]
    starts: tuple[BornStartDiagnostic, ...]
    confirmation_parameter_scaled_max_deltas: tuple[float, float]
    failure_reasons: tuple[str, ...]

    @property
    def observations(self) -> tuple[BornObservationDiagnostic, ...]:
        return self.starts[0].observations


@dataclass(frozen=True, slots=True)
class AqueousMiacRowDiagnostic:
    row_id: str
    salt: str
    molality_mol_per_kg: float
    observed_gamma_pm_m: float
    modeled_log_gamma_pm_m: float
    modeled_gamma_pm_m: float
    raw_error: float
    scaled_residual: float
    local_log_derivative: tuple[float, ...]
    reference_molality_mol_per_kg: float
    reference_convergence_error: float
    reference_derivative_convergence_error: float
    provider_fingerprint: str


@dataclass(frozen=True, slots=True)
class AqueousStageStartDiagnostic:
    stage: str
    name: str
    termination: str
    solution_usable: bool
    initial_cost: float
    final_cost: float
    iterations: int
    parameters: tuple[float, ...]
    rows: tuple[AqueousMiacRowDiagnostic, ...]
    singular_values: tuple[float, ...]
    rank_threshold: float
    rank: int
    condition_number: float
    least_sensitive_direction: tuple[float, ...]
    complete_columns: bool
    active_bounds: tuple[bool, ...]
    solver_converged: bool
    numerically_valid: bool
    failure_reasons: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class FigielStagedCycleDiagnostic:
    cycle_index: int
    born_starts: tuple[BornStartDiagnostic, ...]
    solvation_factor_starts: tuple[AqueousStageStartDiagnostic, ...]
    aqueous_kij_starts: tuple[AqueousStageStartDiagnostic, ...]
    born_diameters_angstrom: tuple[float, ...]
    water_solvation_factor: float
    aqueous_kij: tuple[float, ...]
    scaled_max_delta_from_previous: float | None
    cycle_converged: bool


@dataclass(frozen=True, slots=True)
class FigielStagedAqueousRecoveryResult:
    specification_id: str
    compiled_problem_identities: tuple[tuple[str, ...], ...]
    provider_fingerprints: tuple[str, ...]
    solver_converged: bool
    numerically_converged: bool
    physically_valid: bool
    workflow_valid: bool
    scientifically_valid: bool
    predictive_status: str
    born_diameters_angstrom: tuple[float, ...]
    water_solvation_factor: float
    aqueous_kij: tuple[float, ...]
    published_aqueous_kij: tuple[float, ...]
    maximum_published_kij_difference: float
    pooled_miac_rmse: float
    per_salt_miac_rmse: tuple[tuple[str, float], ...]
    per_salt_miac_max_abs_error: tuple[tuple[str, float], ...]
    first_predicted_miac: tuple[tuple[str, float], ...]
    input_row_ids: tuple[str, ...]
    evaluated_row_ids: tuple[str, ...]
    failed_row_ids: tuple[str, ...]
    cycles: tuple[FigielStagedCycleDiagnostic, ...]
    final_rows: tuple[AqueousMiacRowDiagnostic, ...]
    failure_reasons: tuple[str, ...]


def _born_native_payload(
    specification: BornDiameterTracerSpecification,
    *,
    expected_fingerprints: tuple[str, ...] | None = None,
    starts: tuple[tuple[float, ...], ...] | None = None,
    staged: bool = False,
) -> tuple[object, ...]:
    identity = (
        specification.specification_id,
        specification.source_validation_commit,
        specification.source_validation_tree,
        specification.source_ledger_sha256,
        specification.source_parameter_packet_sha256,
        specification.source_metadata_sha256,
        specification.packaged_targets_sha256,
        specification.source_doi,
        specification.source_si_doi,
        specification.source_locator,
        specification.source_basis,
        "K",
        "Pa",
        "J/mol",
        "angstrom",
        "x-process at infinite dilution",
        PROVIDER_CAPSULE,
        DIAMETER_TRANSFORM,
        BORN_RESIDUAL,
        BORN_JACOBIAN,
        specification.ceres_linear_solver,
        specification.ceres_logging,
        *(('figiel-staged-aqueous-recovery',) if staged else ()),
    )
    fingerprints = expected_fingerprints or tuple(
        target.expected_provider_fingerprint for target in specification.targets
    )
    if len(fingerprints) != len(specification.targets):
        raise ValueError("Born Provider fingerprints must contain five entries")
    targets = tuple(
        (
            target.target_id,
            target.ion_label,
            target.active_component_id,
            target.counterion_component_id,
            target.target_j_per_mol,
            target.published_diameter_angstrom,
            fingerprint,
        )
        for target, fingerprint in zip(
            specification.targets, fingerprints, strict=True
        )
    )
    return (
        identity,
        targets,
        specification.temperature_k,
        specification.pressure_pa,
        specification.reference_molality_mol_per_kg,
        specification.reference_convergence_error_max,
        specification.diameter_origin_angstrom,
        specification.diameter_scale_angstrom,
        specification.diameter_bounds_angstrom,
        specification.scaled_bounds,
        starts or specification.start_diameters_angstrom,
        specification.max_num_iterations,
        specification.function_tolerance,
        specification.gradient_tolerance,
        specification.parameter_tolerance,
        specification.rank_threshold_multiplier,
    )


def _native_payload(
    dataset: PureSaturationDataset,
    specification: PureSaturationFitSpecification,
    provider_fingerprint: str,
) -> tuple[object, ...]:
    identity = (
        dataset.dataset_id,
        dataset.component_id,
        dataset.temperature_unit,
        dataset.pressure_unit,
        dataset.liquid_density_unit,
        dataset.source.source_id,
        dataset.source.citation,
        dataset.source.locator,
        dataset.source.url,
        dataset.source.retrieved_on,
        dataset.source.use_basis,
        dataset.source.transformation,
        dataset.source.units[0][0],
        dataset.source.units[0][1],
        dataset.source.units[1][0],
        dataset.source.units[1][1],
        dataset.source.units[2][0],
        dataset.source.units[2][1],
        dataset.source.data_sha256,
        dataset.source.packaged_data_sha256,
        specification.specification_id,
        *specification.parameter_names,
        *specification.parameter_units,
        *specification.residual_names,
        "mol",
        "m3/mol",
        PROVIDER_CAPSULE,
        provider_fingerprint,
        PARAMETER_TRANSFORM,
        LIQUID_VOLUME_TRANSFORM,
        VAPOR_VOLUME_TRANSFORM,
        REPORTING_PRESSURE_TRANSFORM,
        specification.ceres_linear_solver,
        specification.ceres_logging,
    )
    return (
        identity,
        tuple(_row_payload(row) for row in dataset.training_rows),
        specification.start,
        specification.lower_bounds,
        specification.upper_bounds,
        specification.parameter_scales,
        specification.fixed_amount_mol,
        specification.molar_mass_kg_per_mol,
        specification.residual_weights,
        specification.liquid_volume_bounds_m3,
        specification.vapor_volume_bounds_m3,
        specification.topology_relative_separation_min,
        specification.max_num_iterations,
        specification.function_tolerance,
        specification.gradient_tolerance,
        specification.parameter_tolerance,
        specification.reporting_pressure_bounds_pa,
        specification.confirmation_liquid_volume_start_multiplier,
        specification.confirmation_vapor_volume_start_multiplier,
        specification.confirmation_parameter_scaled_max_delta,
        specification.confirmation_cost_relative_delta,
        specification.reporting_pressure_scaled_residual_max,
        specification.reporting_chemical_potential_residual_max,
        specification.ceres_num_threads,
    )


def _active_bound(value: float, lower: float, upper: float) -> str | None:
    tolerance = 1.0e-8 * (1.0 + max(abs(lower), abs(upper)))
    if abs(value - lower) <= tolerance:
        return "lower"
    if abs(value - upper) <= tolerance:
        return "upper"
    return None


def _reporting_row_diagnostic(
    source: SaturationObservation,
    training_ids: frozenset[str],
    held_out_ids: frozenset[str],
    stress_ids: frozenset[str],
    specification: PureSaturationFitSpecification,
    native_row: tuple[object, ...],
) -> ReportingRowDiagnostic:
    if str(native_row[0]) != source.row_id or str(native_row[1]) != source.source_id:
        raise RuntimeError("native reporting row identity did not match the immutable dataset")
    predicted_pressure = float(native_row[5])
    predicted_density = float(native_row[6])
    liquid_volume = float(native_row[7])
    vapor_volume = float(native_row[8])
    liquid_slope = float(native_row[9])
    vapor_slope = float(native_row[10])
    raw = tuple(float(value) for value in native_row[11])
    termination = str(native_row[12])
    usable = bool(native_row[13])
    native_failure_reason = str(native_row[14]).strip()
    reasons = [native_failure_reason] if native_failure_reason else []
    finite = all(
        math.isfinite(value)
        for value in (
            predicted_pressure,
            predicted_density,
            liquid_volume,
            vapor_volume,
            liquid_slope,
            vapor_slope,
            *raw,
        )
    )
    topology_ok = (
        finite
        and liquid_volume > 0.0
        and vapor_volume > 0.0
        and liquid_volume < vapor_volume
        and (vapor_volume - liquid_volume) / vapor_volume
        > specification.topology_relative_separation_min
    )
    if termination != "CONVERGENCE":
        reasons.append(f"reporting Ceres termination was {termination}")
    if not usable:
        reasons.append("reporting Ceres solution was unusable")
    if not finite:
        reasons.append("reporting diagnostics were nonfinite")
    if not topology_ok:
        reasons.append("reporting phases failed the topology separation gate")
    if finite and (liquid_slope <= 0.0 or vapor_slope <= 0.0):
        reasons.append("reporting phase was mechanically unstable")
    if finite and (
        abs(raw[0]) / source.pressure_pa
        > specification.reporting_pressure_scaled_residual_max
        or abs(raw[1]) / source.pressure_pa
        > specification.reporting_pressure_scaled_residual_max
    ):
        reasons.append("reporting scaled pressure closure exceeded its threshold")
    if finite and abs(raw[2]) > specification.reporting_chemical_potential_residual_max:
        reasons.append("reporting chemical-potential closure exceeded its threshold")

    liquid_molar_density = (
        specification.fixed_amount_mol / liquid_volume
        if math.isfinite(liquid_volume) and liquid_volume > 0.0
        else math.nan
    )
    vapor_molar_density = (
        specification.fixed_amount_mol / vapor_volume
        if math.isfinite(vapor_volume) and vapor_volume > 0.0
        else math.nan
    )
    memberships = (
        source.row_id in training_ids,
        source.row_id in held_out_ids,
        source.row_id in stress_ids,
    )
    if sum(memberships) != 1:
        raise RuntimeError("reporting row did not belong to exactly one immutable partition")
    partition = ("training", "held_out", "stress")[memberships.index(True)]
    return ReportingRowDiagnostic(
        row_id=source.row_id,
        temperature_k=source.temperature_k,
        training=source.row_id in training_ids,
        partition=partition,
        observed_pressure_pa=source.pressure_pa,
        predicted_pressure_pa=predicted_pressure,
        pressure_relative_error=(predicted_pressure - source.pressure_pa) / source.pressure_pa,
        observed_liquid_density_kg_m3=source.liquid_density_kg_m3,
        predicted_liquid_density_kg_m3=predicted_density,
        liquid_density_relative_error=(predicted_density - source.liquid_density_kg_m3)
        / source.liquid_density_kg_m3,
        liquid_volume_m3=liquid_volume,
        vapor_volume_m3=vapor_volume,
        liquid_molar_density_mol_m3=liquid_molar_density,
        vapor_molar_density_mol_m3=vapor_molar_density,
        liquid_mass_density_kg_m3=(
            specification.molar_mass_kg_per_mol * liquid_molar_density
        ),
        vapor_mass_density_kg_m3=(
            specification.molar_mass_kg_per_mol * vapor_molar_density
        ),
        liquid_stability_slope=liquid_slope,
        vapor_stability_slope=vapor_slope,
        raw_equilibrium_residuals=raw,
        termination=termination,
        solution_usable=usable,
        physically_valid=not reasons,
        failure_reasons=tuple(reasons),
    )


def fit_pure_saturation(
    *,
    model: object,
    dataset: PureSaturationDataset,
    specification: PureSaturationFitSpecification,
) -> PureSaturationFitResult:
    if type(dataset) is not PureSaturationDataset:
        raise TypeError("dataset must be an exact PureSaturationDataset")
    if type(specification) is not PureSaturationFitSpecification:
        raise TypeError("specification must be an exact PureSaturationFitSpecification")
    if (
        dataset.component_id != specification.component_id
        or dataset.dataset_id != specification.dataset_id
        or dataset.source.source_id != specification.source_id
        or dataset.training_temperatures_k != specification.training_temperatures_k
    ):
        raise ValueError("dataset and specification identities do not match")
    capsule = native_sdk(model)
    provider_fingerprint = getattr(model, "parameter_fingerprint", None)
    if not isinstance(provider_fingerprint, str) or not provider_fingerprint:
        raise ValueError("model must expose a nonblank provider parameter_fingerprint")
    if provider_fingerprint != specification.expected_provider_fingerprint:
        raise ValueError("model fingerprint does not match the immutable component specification")
    payload = _native_payload(dataset, specification, provider_fingerprint)
    reporting_payload = tuple(_row_payload(row) for row in dataset.rows)
    (
        termination_native,
        solution_usable_native,
        initial_cost_native,
        final_cost_native,
        iterations_native,
        variables_native,
        _residuals_native,
        _jacobian_native,
        training_rows_native,
        full_singular_values_native,
        full_rank_native,
        full_condition_native,
        parameter_singular_values_native,
        parameter_rank_native,
        parameter_condition_native,
        complete_columns_native,
        parameter_delta_native,
        cost_delta_native,
        confirmation_termination_native,
        confirmation_usable_native,
        reporting_rows_native,
        observed_fingerprint_native,
        compiled_identity_native,
        native_failure_reason_native,
    ) = _native.solve(capsule, payload, reporting_payload)
    if tuple(compiled_identity_native) != payload[0]:
        raise RuntimeError("compiled problem identity did not round-trip from the native solve")
    variables = tuple(float(value) for value in variables_native)
    final_parameters = tuple(
        start + scale * transformed
        for start, scale, transformed in zip(
            specification.start, specification.parameter_scales, variables[:3], strict=True
        )
    )
    parameters = tuple(
        ParameterDiagnostic(
            name=name,
            unit=unit,
            start=start,
            final=final,
            movement=final - start,
            lower_bound=lower,
            upper_bound=upper,
            active_bound=_active_bound(final, lower, upper),
        )
        for name, unit, start, final, lower, upper in zip(
            specification.parameter_names,
            specification.parameter_units,
            specification.start,
            final_parameters,
            specification.lower_bounds,
            specification.upper_bounds,
            strict=True,
        )
    )
    native_failure_reason = str(native_failure_reason_native).strip()
    training_rows = tuple(
        TrainingRowDiagnostic(
            row_id=str(native_row[0]),
            temperature_k=float(native_row[2]),
            observed_pressure_pa=source.pressure_pa,
            observed_liquid_density_kg_m3=source.liquid_density_kg_m3,
            liquid_volume_m3=float(native_row[3]),
            vapor_volume_m3=float(native_row[4]),
            liquid_molar_density_mol_m3=specification.fixed_amount_mol / float(native_row[3]),
            vapor_molar_density_mol_m3=specification.fixed_amount_mol / float(native_row[4]),
            liquid_mass_density_kg_m3=specification.molar_mass_kg_per_mol
            / float(native_row[3]),
            vapor_mass_density_kg_m3=specification.molar_mass_kg_per_mol
            / float(native_row[4]),
            liquid_pressure_pa=float(native_row[5]),
            vapor_pressure_pa=float(native_row[6]),
            liquid_chemical_potential_over_rt=float(native_row[7]),
            vapor_chemical_potential_over_rt=float(native_row[8]),
            liquid_stability_slope=float(native_row[9]),
            vapor_stability_slope=float(native_row[10]),
            raw_residuals=tuple(float(value) for value in native_row[11]),
            scaled_residuals=tuple(float(value) for value in native_row[12]),
        )
        for source, native_row in zip(dataset.training_rows, training_rows_native)
        if str(native_row[0]) == source.row_id and str(native_row[1]) == source.source_id
    )
    if training_rows_native and (
        len(training_rows_native) != len(dataset.training_rows)
        or len(training_rows) != len(dataset.training_rows)
    ):
        raise RuntimeError("native training row identity did not match the immutable dataset")
    training_ids = frozenset(dataset.training_row_ids)
    held_out_ids = frozenset(row.row_id for row in dataset.held_out_rows)
    stress_ids = frozenset(row.row_id for row in dataset.stress_rows)
    reporting_rows = [
        _reporting_row_diagnostic(
            source,
            training_ids,
            held_out_ids,
            stress_ids,
            specification,
            native_row,
        )
        for source, native_row in zip(dataset.rows, reporting_rows_native)
    ]
    if reporting_rows_native and len(reporting_rows_native) != len(dataset.rows):
        raise RuntimeError("native reporting row identity did not match the immutable dataset")
    jacobian = JacobianDiagnostics(
        complete_columns=bool(complete_columns_native),
        full_singular_values=tuple(float(value) for value in full_singular_values_native),
        full_rank=int(full_rank_native),
        full_condition_number=float(full_condition_native),
        parameter_singular_values=tuple(
            float(value) for value in parameter_singular_values_native
        ),
        parameter_rank=int(parameter_rank_native),
        parameter_condition_number=float(parameter_condition_native),
    )
    termination = str(termination_native)
    usable = bool(solution_usable_native)
    initial_cost = float(initial_cost_native)
    final_cost = float(final_cost_native)
    bounds_respected = all(
        item.lower_bound <= item.final <= item.upper_bound for item in parameters
    )
    fitted_parameter_count = len(parameters)
    parameter_columns_full_rank = jacobian.parameter_rank == fitted_parameter_count
    solver_converged = (
        termination == "CONVERGENCE"
        and usable
        and math.isfinite(initial_cost)
        and math.isfinite(final_cost)
        and final_cost <= initial_cost
        and jacobian.complete_columns
        and parameter_columns_full_rank
        and bounds_respected
        and not native_failure_reason
    )
    confirmation_termination = str(confirmation_termination_native)
    confirmation_usable = bool(confirmation_usable_native)
    parameter_delta = float(parameter_delta_native)
    cost_delta = float(cost_delta_native)
    numerically_converged = (
        solver_converged
        and confirmation_termination == "CONVERGENCE"
        and confirmation_usable
        and parameter_delta <= specification.confirmation_parameter_scaled_max_delta
        and cost_delta <= specification.confirmation_cost_relative_delta
    )
    reporting_tuple = tuple(reporting_rows)
    acceptance_reporting_rows = tuple(
        row for row in reporting_tuple if row.partition != "stress"
    )
    physical_valid = (
        solver_converged
        and all(row.liquid_volume_m3 < row.vapor_volume_m3 for row in training_rows)
        and all(row.liquid_stability_slope > 0.0 for row in training_rows)
        and all(row.vapor_stability_slope > 0.0 for row in training_rows)
        and all(row.physically_valid for row in acceptance_reporting_rows)
    )
    failure_reasons: list[str] = []
    if native_failure_reason:
        failure_reasons.append(native_failure_reason)
    failure_reasons.extend(
        f"{row.row_id}: {reason}"
        for row in acceptance_reporting_rows
        for reason in row.failure_reasons
    )
    if not parameter_columns_full_rank:
        failure_reasons.append(
            "training parameter Jacobian is rank deficient: "
            f"{jacobian.parameter_rank} of {fitted_parameter_count} fitted parameter columns"
        )
    if not solver_converged:
        failure_reasons.append("training solver convergence gate failed")
    if not numerically_converged:
        failure_reasons.append("confirmation solve numerical convergence gate failed")
    if not physical_valid:
        failure_reasons.append("training or reporting physical validity gate failed")
    observed_fingerprint = str(observed_fingerprint_native)
    if observed_fingerprint and observed_fingerprint != getattr(
        model, "parameter_fingerprint", None
    ):
        failure_reasons.append("provider source fingerprint did not match the supplied model")
        physical_valid = False
    return PureSaturationFitResult(
        component_id=dataset.component_id,
        dataset_id=dataset.dataset_id,
        specification_id=specification.specification_id,
        provider_fingerprint=observed_fingerprint,
        compiled_problem_identity=tuple(str(value) for value in compiled_identity_native),
        solver_converged=solver_converged,
        numerically_converged=numerically_converged,
        physically_valid=physical_valid,
        predictive_status=PREDICTIVE_STATUS,
        termination=termination,
        solution_usable=usable,
        initial_cost=initial_cost,
        final_cost=final_cost,
        iterations=int(iterations_native),
        parameters=parameters,
        jacobian=jacobian,
        training_rows=training_rows,
        reporting_rows=reporting_tuple,
        confirmation_termination=confirmation_termination,
        confirmation_solution_usable=confirmation_usable,
        confirmation_parameter_scaled_max_delta=parameter_delta,
        confirmation_cost_relative_delta=cost_delta,
        failure_reasons=tuple(failure_reasons),
    )


def _born_observations(
    native_rows: tuple[object, ...],
    specification: BornDiameterTracerSpecification,
) -> tuple[BornObservationDiagnostic, ...]:
    if len(native_rows) != len(specification.targets):
        raise RuntimeError("native Born observations did not match the five-target contract")
    return tuple(
        BornObservationDiagnostic(
            target_id=target.target_id,
            ion_label=target.ion_label,
            target_j_per_mol=target.target_j_per_mol,
            modeled_j_per_mol=float(row[0]),
            derivative_j_per_mol_per_angstrom=float(row[1]),
            raw_error_j_per_mol=float(row[2]),
            scaled_residual=float(row[3]),
            scaled_jacobian=float(row[4]),
            reference_molality_mol_per_kg=float(row[5]),
            reference_convergence_error=float(row[6]),
            provider_fingerprint=str(row[7]),
        )
        for target, row in zip(specification.targets, native_rows, strict=True)
    )


def _born_start_diagnostic(
    native_start: tuple[object, ...],
    specification: BornDiameterTracerSpecification,
) -> BornStartDiagnostic:
    (
        name_native,
        termination_native,
        usable_native,
        initial_cost_native,
        final_cost_native,
        iterations_native,
        transformed_native,
        residuals_native,
        jacobian_native,
        rows_native,
        singular_values_native,
        rank_native,
        condition_native,
        complete_columns_native,
        failure_native,
    ) = native_start
    transformed = tuple(float(value) for value in transformed_native)
    residuals = tuple(float(value) for value in residuals_native)
    jacobian = tuple(float(value) for value in jacobian_native)
    singular_values = tuple(float(value) for value in singular_values_native)
    if not (
        len(transformed) == len(residuals) == len(singular_values) == 5
        and len(jacobian) == 25
    ):
        raise RuntimeError("native Born result dimensions did not match the 5 x 5 contract")
    final_diameters = tuple(
        specification.diameter_origin_angstrom
        + specification.diameter_scale_angstrom * value
        for value in transformed
    )
    active_tolerance = math.sqrt(math.ulp(1.0)) * max(
        1.0,
        abs(specification.scaled_bounds[0]),
        abs(specification.scaled_bounds[1]),
    )
    inactive_bounds = all(
        min(
            value - specification.scaled_bounds[0],
            specification.scaled_bounds[1] - value,
        )
        > active_tolerance
        for value in transformed
    )
    observations = _born_observations(tuple(rows_native), specification)
    initial_cost = float(initial_cost_native)
    final_cost = float(final_cost_native)
    condition_number = float(condition_native)
    rank = int(rank_native)
    complete_columns = bool(complete_columns_native)
    native_failure = str(failure_native).strip()
    finite = all(
        math.isfinite(value)
        for value in (
            initial_cost,
            final_cost,
            condition_number,
            *transformed,
            *residuals,
            *jacobian,
            *singular_values,
            *(row.modeled_j_per_mol for row in observations),
            *(row.derivative_j_per_mol_per_angstrom for row in observations),
        )
    )
    reasons: list[str] = []
    if str(termination_native) != "CONVERGENCE":
        reasons.append(f"Ceres termination was {termination_native}")
    if not bool(usable_native):
        reasons.append("Ceres solution was unusable")
    if not finite:
        reasons.append("Born solution or Jacobian was nonfinite")
    if final_cost > initial_cost + math.ulp(max(1.0, abs(initial_cost))):
        reasons.append("Born solve increased cost")
    if not complete_columns:
        reasons.append("Born Jacobian columns were incomplete")
    if rank != 5:
        reasons.append(f"Born Jacobian rank was {rank} of 5")
    if not inactive_bounds:
        reasons.append("Born solution had an active or violated bound")
    if native_failure:
        reasons.append(native_failure)
    rank_threshold = (
        singular_values[0]
        * 5.0
        * math.ulp(1.0)
        * specification.rank_threshold_multiplier
    )
    return BornStartDiagnostic(
        name=str(name_native),
        termination=str(termination_native),
        solution_usable=bool(usable_native),
        initial_cost=initial_cost,
        final_cost=final_cost,
        iterations=int(iterations_native),
        transformed_parameters=transformed,
        final_diameters_angstrom=final_diameters,
        observations=observations,
        singular_values=singular_values,
        rank_threshold=rank_threshold,
        rank=rank,
        condition_number=condition_number,
        complete_columns=complete_columns,
        inactive_bounds=inactive_bounds,
        solver_converged=not reasons,
        failure_reasons=tuple(reasons),
    )


def fit_figiel_born_diameters(*, models: tuple[object, ...]) -> BornDiameterFitResult:
    specification = FIGIEL_BORN_DIAMETER_TRACER_V1
    if type(models) is not tuple or len(models) != 5:
        raise TypeError("models must be the exact ordered five-model tuple")
    for model, target in zip(models, specification.targets, strict=True):
        component_ids = getattr(model, "component_ids", None)
        if tuple(component_ids or ()) != target.component_order:
            raise ValueError("model component order does not match the immutable Born target")
        fingerprint = getattr(model, "parameter_fingerprint", None)
        if fingerprint != target.expected_provider_fingerprint:
            raise ValueError("model fingerprint does not match the immutable Born target")
    capsules = tuple(native_sdk(model) for model in models)
    payload = _born_native_payload(specification)
    starts_native, compiled_identity_native = _native.solve_born(capsules, payload)
    if tuple(compiled_identity_native) != payload[0]:
        raise RuntimeError("compiled Born problem identity did not round-trip")
    starts = tuple(
        _born_start_diagnostic(tuple(native_start), specification)
        for native_start in starts_native
    )
    if tuple(start.name for start in starts) != ("primary", "lower", "upper"):
        raise RuntimeError("native Born starts did not match the frozen schedule")
    primary = starts[0]
    parameters = tuple(
        BornParameterDiagnostic(
            ion_label=target.ion_label,
            active_component_id=target.active_component_id,
            final_diameter_angstrom=diameter,
            published_diameter_angstrom=target.published_diameter_angstrom,
            published_delta_angstrom=diameter - target.published_diameter_angstrom,
            lower_bound_angstrom=specification.diameter_bounds_angstrom[0],
            upper_bound_angstrom=specification.diameter_bounds_angstrom[1],
            scaled_lower_bound_distance=transformed - specification.scaled_bounds[0],
            scaled_upper_bound_distance=specification.scaled_bounds[1] - transformed,
            active_bound=not (
                min(
                    transformed - specification.scaled_bounds[0],
                    specification.scaled_bounds[1] - transformed,
                )
                > math.sqrt(math.ulp(1.0))
                * max(1.0, *(abs(value) for value in specification.scaled_bounds))
            ),
        )
        for target, diameter, transformed in zip(
            specification.targets,
            primary.final_diameters_angstrom,
            primary.transformed_parameters,
            strict=True,
        )
    )
    confirmation_deltas = tuple(
        max(
            abs(value - reference)
            for value, reference in zip(
                confirmation.transformed_parameters,
                primary.transformed_parameters,
                strict=True,
            )
        )
        for confirmation in starts[1:]
    )
    solver_converged = primary.solver_converged
    numerical_converged = (
        solver_converged
        and all(start.solver_converged for start in starts)
        and all(
            max(abs(row.scaled_residual) for row in start.observations)
            <= specification.scaled_residual_max
            for start in starts
        )
        and all(
            delta <= specification.confirmation_parameter_scaled_max_delta
            for delta in confirmation_deltas
        )
    )
    expected_fingerprints = tuple(
        target.expected_provider_fingerprint for target in specification.targets
    )
    observed_fingerprints = tuple(row.provider_fingerprint for row in primary.observations)
    workflow_valid = (
        numerical_converged
        and observed_fingerprints == expected_fingerprints
        and all(
            row.reference_molality_mol_per_kg
            == specification.reference_molality_mol_per_kg
            and row.reference_convergence_error
            <= specification.reference_convergence_error_max
            for start in starts
            for row in start.observations
        )
    )
    scientifically_valid = workflow_valid and all(
        abs(row.raw_error_j_per_mol)
        <= specification.observable_round_trip_j_per_mol
        for start in starts
        for row in start.observations
    )
    failure_reasons = [
        f"{start.name}: {reason}"
        for start in starts
        for reason in start.failure_reasons
    ]
    if not numerical_converged:
        failure_reasons.append("three-start numerical confirmation gate failed")
    if not workflow_valid:
        failure_reasons.append("source-bound workflow identity or reference gate failed")
    if not scientifically_valid:
        failure_reasons.append("source-observable reproduction gate failed")
    return BornDiameterFitResult(
        specification_id=specification.specification_id,
        compiled_problem_identity=tuple(str(value) for value in compiled_identity_native),
        provider_fingerprints=observed_fingerprints,
        solver_converged=solver_converged,
        numerically_converged=numerical_converged,
        workflow_valid=workflow_valid,
        scientifically_valid=scientifically_valid,
        predictive_status=PREDICTIVE_STATUS,
        parameters=parameters,
        starts=starts,
        confirmation_parameter_scaled_max_deltas=(
            float(confirmation_deltas[0]),
            float(confirmation_deltas[1]),
        ),
        failure_reasons=tuple(failure_reasons),
    )


def _figiel_trial_bundle(
    born_diameters_angstrom: tuple[float, ...],
    water_solvation_factor: float,
    aqueous_kij: tuple[float, ...],
) -> ParameterBundle:
    if len(born_diameters_angstrom) != 5 or len(aqueous_kij) != 11:
        raise ValueError("staged Figiel parameter dimensions must be 5 + 1 + 11")
    catalog = ParameterBundle.from_catalog(
        "figiel-2025-reference-electrolytes", version=1
    )
    born_by_component = {
        target.active_component_id: value
        for target, value in zip(
            FIGIEL_BORN_DIAMETER_TRACER_V1.targets,
            born_diameters_angstrom,
            strict=True,
        )
    }
    kij_by_pair = {
        frozenset(pair): value
        for pair, value in zip(
            FIGIEL_AQUEOUS_KIJ_COORDINATES, aqueous_kij, strict=True
        )
    }
    records = []
    for record in catalog.records:
        if (
            isinstance(record, SingleParameterRecord)
            and record.family == "born_diameter"
            and record.component_id in born_by_component
        ):
            records.append(
                replace(
                    record,
                    value=(
                        born_by_component[record.component_id]
                        * unit_registry.angstrom
                    ),
                )
            )
        elif (
            isinstance(record, SingleParameterRecord)
            and record.record_id == "water-solvation-factor"
        ):
            records.append(replace(record, value=water_solvation_factor))
        elif isinstance(record, PairParameterRecord) and frozenset(
            (record.component_id_a, record.component_id_b)
        ) in kij_by_pair:
            records.append(
                replace(
                    record,
                    value=kij_by_pair[
                        frozenset((record.component_id_a, record.component_id_b))
                    ],
                )
            )
        else:
            records.append(record)
    return ParameterBundle.from_records(
        bundle_id="figiel-staged-aqueous-trial",
        bundle_version=1,
        purpose="user-provided",
        sources=catalog.sources,
        domains=catalog.domains,
        components=catalog.components,
        singles=(record for record in records if isinstance(record, SingleParameterRecord)),
        pairs=(record for record in records if isinstance(record, PairParameterRecord)),
        sites=(record for record in records if isinstance(record, SiteRecord)),
        associations=(
            record
            for record in records
            if isinstance(record, AssociationParameterRecord)
        ),
        correlations=(
            record
            for record in records
            if isinstance(
                record,
                (ConstantCorrelation, ConstantPlusSumOfExponentialsCorrelation),
            )
        ),
        models=(record for record in records if isinstance(record, ModelParameterRecord)),
    )


def _figiel_models(
    born_diameters_angstrom: tuple[float, ...],
    water_solvation_factor: float,
    aqueous_kij: tuple[float, ...],
) -> tuple[tuple[EPCSAFT, ...], tuple[EPCSAFT, ...]]:
    bundle = _figiel_trial_bundle(
        born_diameters_angstrom, water_solvation_factor, aqueous_kij
    )
    born_models = tuple(
        EPCSAFT(bundle.select(target.component_order))
        for target in FIGIEL_BORN_DIAMETER_TRACER_V1.targets
    )
    salt_models = tuple(
        EPCSAFT(
            bundle.select(
                (
                    "water",
                    FIGIEL_AQUEOUS_COMPONENTS[salt][0],
                    FIGIEL_AQUEOUS_COMPONENTS[salt][1],
                )
            )
        )
        for salt in FIGIEL_AQUEOUS_SALTS
    )
    return born_models, salt_models


def _aqueous_native_payload(
    specification: FigielStagedAqueousRecoverySpecification,
    *,
    stage: str,
    expected_fingerprints: tuple[str, ...],
    starts: tuple[tuple[str, tuple[float, ...]], ...],
) -> tuple[object, ...]:
    if stage == "solvation_factor":
        observations = specification.stage_b_observations
        parameter_count = 1
        bounds = specification.solvent_factor_bounds
    elif stage == "aqueous_kij":
        observations = specification.observations
        parameter_count = 11
        bounds = specification.kij_bounds
    else:
        raise ValueError("unsupported staged aqueous family")
    rows = tuple(
        (
            row.row_id,
            row.salt,
            FIGIEL_AQUEOUS_SALTS.index(row.salt),
            row.molality_mol_per_kg,
            row.gamma_pm_m,
            (0,) if stage == "solvation_factor" else AQUEOUS_KIJ_COLUMNS[row.salt],
        )
        for row in observations
    )
    identity = (
        specification.specification_id,
        specification.source_validation_commit,
        specification.source_validation_tree,
        specification.source_ledger_sha256,
        specification.source_parameter_packet_sha256,
        specification.source_metadata_sha256,
        specification.source_si_extraction_sha256,
        specification.source_csv_sha256,
        stage,
        "298.15 K",
        "100000 Pa",
        "mol/kg",
        "dimensionless molality-scale mean ionic activity coefficient",
        AQUEOUS_MIAC_RESIDUAL,
        AQUEOUS_MIAC_JACOBIAN,
        "equal row weights; observations are not uncertainties",
        PROVIDER_CAPSULE,
        "DENSE_QR",
        "per-aqueous-start wall-time maximum: 180 s",
    )
    return (
        identity,
        stage,
        rows,
        expected_fingerprints,
        parameter_count,
        bounds,
        starts,
        specification.temperature_k,
        specification.pressure_pa,
        specification.max_num_iterations,
        specification.aqueous_start_wall_time_max_seconds,
        specification.function_tolerance,
        specification.gradient_tolerance,
        specification.parameter_tolerance,
        specification.rank_threshold_multiplier,
    )


def _aqueous_rows(native_rows: tuple[object, ...]) -> tuple[AqueousMiacRowDiagnostic, ...]:
    return tuple(
        AqueousMiacRowDiagnostic(
            row_id=str(row[0]),
            salt=str(row[1]),
            molality_mol_per_kg=float(row[2]),
            observed_gamma_pm_m=float(row[3]),
            modeled_log_gamma_pm_m=float(row[4]),
            modeled_gamma_pm_m=float(row[5]),
            raw_error=float(row[5]) - float(row[3]),
            scaled_residual=float(row[6]),
            local_log_derivative=tuple(float(value) for value in row[7]),
            reference_molality_mol_per_kg=float(row[8]),
            reference_convergence_error=float(row[9]),
            reference_derivative_convergence_error=float(row[10]),
            provider_fingerprint=str(row[11]),
        )
        for row in native_rows
    )


def _aqueous_start_diagnostic(
    native_start: tuple[object, ...],
    *,
    stage: str,
    bounds: tuple[float, float],
    expected_rank: int,
) -> AqueousStageStartDiagnostic:
    (
        name_native,
        termination_native,
        usable_native,
        initial_cost_native,
        final_cost_native,
        iterations_native,
        parameters_native,
        residuals_native,
        jacobian_native,
        rows_native,
        singular_native,
        rank_threshold_native,
        rank_native,
        condition_native,
        least_native,
        complete_native,
        failure_native,
    ) = native_start
    parameters = tuple(float(value) for value in parameters_native)
    residuals = tuple(float(value) for value in residuals_native)
    jacobian = tuple(float(value) for value in jacobian_native)
    singular_values = tuple(float(value) for value in singular_native)
    least_sensitive = tuple(float(value) for value in least_native)
    rows = _aqueous_rows(tuple(rows_native))
    parameter_count = 1 if stage == "solvation_factor" else 11
    row_count = 21 if stage == "solvation_factor" else 164
    if not (
        len(parameters) == len(singular_values) == len(least_sensitive) == parameter_count
        and len(residuals) == len(rows) == row_count
        and len(jacobian) == row_count * parameter_count
    ):
        raise RuntimeError("native staged aqueous result dimensions did not round-trip")
    active_tolerance = math.sqrt(math.ulp(1.0)) * max(
        1.0, abs(bounds[0]), abs(bounds[1])
    )
    active_bounds = tuple(
        min(value - bounds[0], bounds[1] - value) <= active_tolerance
        for value in parameters
    )
    initial_cost = float(initial_cost_native)
    final_cost = float(final_cost_native)
    condition_number = float(condition_native)
    rank_threshold = float(rank_threshold_native)
    rank = int(rank_native)
    complete_columns = bool(complete_native)
    native_failure = str(failure_native).strip()
    finite = all(
        math.isfinite(value)
        for value in (
            initial_cost,
            final_cost,
            condition_number,
            rank_threshold,
            *parameters,
            *residuals,
            *jacobian,
            *singular_values,
            *least_sensitive,
            *(row.modeled_gamma_pm_m for row in rows),
            *(row.reference_molality_mol_per_kg for row in rows),
            *(row.reference_convergence_error for row in rows),
            *(row.reference_derivative_convergence_error for row in rows),
        )
    )
    termination = str(termination_native)
    solution_usable = bool(usable_native)
    solver_converged = termination == "CONVERGENCE" and solution_usable
    reasons: list[str] = []
    if termination != "CONVERGENCE":
        reasons.append(f"Ceres termination was {termination}")
    if not solution_usable:
        reasons.append("Ceres solution was unusable")
    if not finite:
        reasons.append("stage solution or exact Jacobian was nonfinite")
    if final_cost > initial_cost + math.ulp(max(1.0, abs(initial_cost))):
        reasons.append("stage solve increased cost")
    if not complete_columns:
        reasons.append("stage Jacobian columns were incomplete")
    if rank != expected_rank:
        reasons.append(f"stage Jacobian rank was {rank} of {expected_rank}")
    if native_failure:
        reasons.append(native_failure)
    return AqueousStageStartDiagnostic(
        stage=stage,
        name=str(name_native),
        termination=termination,
        solution_usable=solution_usable,
        initial_cost=initial_cost,
        final_cost=final_cost,
        iterations=int(iterations_native),
        parameters=parameters,
        rows=rows,
        singular_values=singular_values,
        rank_threshold=rank_threshold,
        rank=rank,
        condition_number=condition_number,
        least_sensitive_direction=least_sensitive,
        complete_columns=complete_columns,
        active_bounds=active_bounds,
        solver_converged=solver_converged,
        numerically_valid=solver_converged and not reasons,
        failure_reasons=tuple(reasons),
    )


def _solve_aqueous_stage(
    salt_models: tuple[EPCSAFT, ...],
    specification: FigielStagedAqueousRecoverySpecification,
    *,
    stage: str,
    starts: tuple[tuple[str, tuple[float, ...]], ...],
) -> tuple[tuple[AqueousStageStartDiagnostic, ...], tuple[str, ...]]:
    fingerprints = tuple(model.parameter_fingerprint for model in salt_models)
    payload = _aqueous_native_payload(
        specification,
        stage=stage,
        expected_fingerprints=fingerprints,
        starts=starts,
    )
    native_starts, compiled_identity = _native.solve_figiel_aqueous(
        tuple(native_sdk(model) for model in salt_models), payload
    )
    if tuple(compiled_identity) != payload[0]:
        raise RuntimeError("compiled staged aqueous identity did not round-trip")
    bounds = (
        specification.solvent_factor_bounds
        if stage == "solvation_factor"
        else specification.kij_bounds
    )
    expected_rank = 1 if stage == "solvation_factor" else 11
    return (
        tuple(
            _aqueous_start_diagnostic(
                tuple(native_start),
                stage=stage,
                bounds=bounds,
                expected_rank=expected_rank,
            )
            for native_start in native_starts
        ),
        tuple(str(value) for value in compiled_identity),
    )


def _solve_staged_born(
    born_models: tuple[EPCSAFT, ...],
    *,
    starts: tuple[tuple[float, ...], ...],
) -> tuple[tuple[BornStartDiagnostic, ...], tuple[str, ...]]:
    specification = FIGIEL_BORN_DIAMETER_TRACER_V1
    fingerprints = tuple(model.parameter_fingerprint for model in born_models)
    payload = _born_native_payload(
        specification,
        expected_fingerprints=fingerprints,
        starts=starts,
        staged=True,
    )
    native_starts, compiled_identity = _native.solve_born(
        tuple(native_sdk(model) for model in born_models), payload
    )
    if tuple(compiled_identity) != payload[0]:
        raise RuntimeError("compiled staged Born identity did not round-trip")
    return (
        tuple(
            _born_start_diagnostic(tuple(native_start), specification)
            for native_start in native_starts
        ),
        tuple(str(value) for value in compiled_identity),
    )


def _start_max_delta(starts: tuple[object, ...], attribute: str) -> float:
    reference = tuple(getattr(starts[0], attribute))
    return max(
        (
            abs(value - expected)
            for start in starts[1:]
            for value, expected in zip(
                tuple(getattr(start, attribute)), reference, strict=True
            )
        ),
        default=0.0,
    )


def fit_figiel_staged_aqueous_parameters() -> FigielStagedAqueousRecoveryResult:
    specification = FIGIEL_STAGED_AQUEOUS_RECOVERY_V1
    born = tuple(
        target.published_diameter_angstrom
        for target in FIGIEL_BORN_DIAMETER_TRACER_V1.targets
    )
    solvation_factor = 1.5
    aqueous_kij = FIGIEL_AQUEOUS_PUBLISHED_KIJ
    cycles: list[FigielStagedCycleDiagnostic] = []
    identities: list[tuple[str, ...]] = []

    for cycle_index in range(specification.max_confirmation_cycles + 1):
        previous = (*born, solvation_factor, *aqueous_kij)
        born_models, _ = _figiel_models(born, solvation_factor, aqueous_kij)
        born_starts_input = (
            FIGIEL_BORN_DIAMETER_TRACER_V1.start_diameters_angstrom
            if cycle_index == 0
            else (born, born, born)
        )
        born_starts, born_identity = _solve_staged_born(
            born_models, starts=born_starts_input
        )
        born = born_starts[0].final_diameters_angstrom

        _, salt_models = _figiel_models(born, solvation_factor, aqueous_kij)
        solvation_starts_input = (
            tuple(
                (name, (value,))
                for name, value in zip(
                    ("primary", "upper"),
                    specification.solvent_factor_starts,
                    strict=True,
                )
            )
            if cycle_index == 0
            else (("cycle", (solvation_factor,)),)
        )
        solvation_starts, solvation_identity = _solve_aqueous_stage(
            salt_models,
            specification,
            stage="solvation_factor",
            starts=solvation_starts_input,
        )
        solvation_factor = solvation_starts[0].parameters[0]

        _, salt_models = _figiel_models(born, solvation_factor, aqueous_kij)
        kij_starts_input = (
            tuple(
                (name, values)
                for name, values in zip(
                    ("primary", "lower", "upper"),
                    specification.kij_starts,
                    strict=True,
                )
            )
            if cycle_index == 0
            else (("cycle", aqueous_kij),)
        )
        kij_starts, kij_identity = _solve_aqueous_stage(
            salt_models,
            specification,
            stage="aqueous_kij",
            starts=kij_starts_input,
        )
        aqueous_kij = kij_starts[0].parameters
        current = (*born, solvation_factor, *aqueous_kij)
        delta = (
            None
            if cycle_index == 0
            else max(abs(value - prior) for value, prior in zip(current, previous, strict=True))
        )
        cycle_converged = (
            delta is not None and delta <= specification.cycle_scaled_max_delta
        )
        cycles.append(
            FigielStagedCycleDiagnostic(
                cycle_index=cycle_index,
                born_starts=born_starts,
                solvation_factor_starts=solvation_starts,
                aqueous_kij_starts=kij_starts,
                born_diameters_angstrom=born,
                water_solvation_factor=solvation_factor,
                aqueous_kij=aqueous_kij,
                scaled_max_delta_from_previous=delta,
                cycle_converged=cycle_converged,
            )
        )
        identities.extend((born_identity, solvation_identity, kij_identity))
        if cycle_converged:
            break

    final_cycle = cycles[-1]
    final_rows = final_cycle.aqueous_kij_starts[0].rows
    input_row_ids = tuple(row.row_id for row in specification.observations)
    evaluated_row_ids = tuple(row.row_id for row in final_rows)
    failed_row_ids = tuple(row_id for row_id in input_row_ids if row_id not in evaluated_row_ids)
    all_born_starts = tuple(start for cycle in cycles for start in cycle.born_starts)
    all_solvation_starts = tuple(
        start for cycle in cycles for start in cycle.solvation_factor_starts
    )
    all_kij_starts = tuple(
        start for cycle in cycles for start in cycle.aqueous_kij_starts
    )
    solver_converged = all(
        start.termination == "CONVERGENCE" and start.solution_usable
        for start in all_born_starts
    ) and all(
        start.solver_converged
        for start in (*all_solvation_starts, *all_kij_starts)
    )
    initial = cycles[0]
    initial_start_agreement = (
        _start_max_delta(initial.born_starts, "transformed_parameters")
        <= specification.cycle_scaled_max_delta
        and _start_max_delta(initial.solvation_factor_starts, "parameters")
        <= specification.cycle_scaled_max_delta
        and _start_max_delta(initial.aqueous_kij_starts, "parameters")
        <= specification.cycle_scaled_max_delta
    )
    numerically_converged = (
        solver_converged
        and all(not start.failure_reasons for start in all_born_starts)
        and all(start.numerically_valid for start in all_solvation_starts)
        and all(start.numerically_valid for start in all_kij_starts)
        and initial_start_agreement
        and final_cycle.cycle_converged
    )
    born_rows = final_cycle.born_starts[0].observations
    physically_valid = all(
        math.isfinite(row.modeled_j_per_mol)
        and row.reference_molality_mol_per_kg
        == FIGIEL_BORN_DIAMETER_TRACER_V1.reference_molality_mol_per_kg
        and row.reference_convergence_error
        <= FIGIEL_BORN_DIAMETER_TRACER_V1.reference_convergence_error_max
        for row in born_rows
    ) and all(
        math.isfinite(row.modeled_gamma_pm_m)
        and row.modeled_gamma_pm_m > 0.0
        and math.isfinite(row.reference_molality_mol_per_kg)
        and math.isfinite(row.reference_convergence_error)
        and math.isfinite(row.reference_derivative_convergence_error)
        and row.provider_fingerprint.startswith("sha256:")
        for start in (*all_solvation_starts, *all_kij_starts)
        for row in start.rows
    )
    final_salt_fingerprints = tuple(
        next(row.provider_fingerprint for row in final_rows if row.salt == salt)
        for salt in FIGIEL_AQUEOUS_SALTS
    )
    final_born_fingerprints = tuple(row.provider_fingerprint for row in born_rows)
    provider_fingerprints = (*final_born_fingerprints, *final_salt_fingerprints)
    workflow_valid = (
        input_row_ids == evaluated_row_ids
        and not failed_row_ids
        and len(final_rows) == 164
        and all(len(cycle.solvation_factor_starts[0].rows) == 21 for cycle in cycles)
        and all(len(cycle.aqueous_kij_starts[0].rows) == 164 for cycle in cycles)
        and all(fingerprint.startswith("sha256:") for fingerprint in provider_fingerprints)
    )

    pooled_rmse = math.sqrt(
        sum(row.raw_error * row.raw_error for row in final_rows) / len(final_rows)
    )
    per_salt_rmse = tuple(
        (
            salt,
            math.sqrt(
                sum(row.raw_error * row.raw_error for row in final_rows if row.salt == salt)
                / sum(row.salt == salt for row in final_rows)
            ),
        )
        for salt in FIGIEL_AQUEOUS_SALTS
    )
    per_salt_max = tuple(
        (
            salt,
            max(abs(row.raw_error) for row in final_rows if row.salt == salt),
        )
        for salt in FIGIEL_AQUEOUS_SALTS
    )
    first_predicted = tuple(
        (
            salt,
            next(row.modeled_gamma_pm_m for row in final_rows if row.salt == salt),
        )
        for salt in FIGIEL_AQUEOUS_SALTS
    )
    maximum_published_difference = max(
        abs(value - published)
        for value, published in zip(
            aqueous_kij, specification.published_kij, strict=True
        )
    )
    born_observable_gate = all(
        abs(row.raw_error_j_per_mol)
        <= FIGIEL_BORN_DIAMETER_TRACER_V1.observable_round_trip_j_per_mol
        for row in born_rows
    )
    observable_gates = (
        pooled_rmse <= specification.pooled_miac_rmse_max
        and all(value <= specification.per_salt_miac_rmse_max for _, value in per_salt_rmse)
        and all(
            value <= specification.per_salt_miac_max_abs_error
            for _, value in per_salt_max
        )
        and all(value < specification.first_predicted_miac_max for _, value in first_predicted)
    )
    scientific_gates = (
        maximum_published_difference <= specification.parameter_comparison_max_abs
        and born_observable_gate
        and observable_gates
    )
    scientifically_valid = (
        solver_converged
        and numerically_converged
        and physically_valid
        and workflow_valid
        and scientific_gates
    )
    failure_reasons = [
        f"cycle {cycle.cycle_index} Born {start.name}: {reason}"
        for cycle in cycles
        for start in cycle.born_starts
        for reason in start.failure_reasons
    ]
    failure_reasons.extend(
        f"cycle {cycle.cycle_index} {start.stage} {start.name}: {reason}"
        for cycle in cycles
        for start in (*cycle.solvation_factor_starts, *cycle.aqueous_kij_starts)
        for reason in start.failure_reasons
    )
    if not initial_start_agreement:
        failure_reasons.append("declared-start agreement gate failed")
    if not final_cycle.cycle_converged:
        failure_reasons.append("three-cycle numerical confirmation gate failed")
    if not physically_valid:
        failure_reasons.append("Provider state or reference diagnostic gate failed")
    if not workflow_valid:
        failure_reasons.append("source-bound workflow identity or row accounting failed")
    if not scientific_gates:
        failure_reasons.append(
            "SOURCE_DESCRIBED_STAGED_RECOVERY_DID_NOT_REPRODUCE_PRINTED_TUPLE"
        )
    return FigielStagedAqueousRecoveryResult(
        specification_id=specification.specification_id,
        compiled_problem_identities=tuple(identities),
        provider_fingerprints=provider_fingerprints,
        solver_converged=solver_converged,
        numerically_converged=numerically_converged,
        physically_valid=physically_valid,
        workflow_valid=workflow_valid,
        scientifically_valid=scientifically_valid,
        predictive_status=PREDICTIVE_STATUS,
        born_diameters_angstrom=born,
        water_solvation_factor=solvation_factor,
        aqueous_kij=aqueous_kij,
        published_aqueous_kij=specification.published_kij,
        maximum_published_kij_difference=maximum_published_difference,
        pooled_miac_rmse=pooled_rmse,
        per_salt_miac_rmse=per_salt_rmse,
        per_salt_miac_max_abs_error=per_salt_max,
        first_predicted_miac=first_predicted,
        input_row_ids=input_row_ids,
        evaluated_row_ids=evaluated_row_ids,
        failed_row_ids=failed_row_ids,
        cycles=tuple(cycles),
        final_rows=final_rows,
        failure_reasons=tuple(failure_reasons),
    )
