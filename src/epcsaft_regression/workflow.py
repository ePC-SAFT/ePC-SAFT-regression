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
    FIGIEL_AQUEOUS_KIJ_V1,
    FIGIEL_BORN_DIAMETER_TRACER_V1,
    FIGIEL_WATER_SOLVATION_FACTOR_V1,
    BornDiameterTracerSpecification,
    FigielAqueousKijSpecification,
    FigielWaterSolvationFactorSpecification,
    PureSaturationDataset,
    PureSaturationFitSpecification,
    SaturationObservation,
)


PROVIDER_CAPSULE = "epcsaft.native_sdk.v1"
PARAMETER_TRANSFORM = "p_j = start_j + parameter_scale_j * z_j"
LIQUID_VOLUME_TRANSFORM = (
    "V_liquid = (molar_mass / observed_liquid_density) * exp(u_liquid)"
)
VAPOR_VOLUME_TRANSFORM = "V_vapor = (R * T / observed_pressure) * exp(u_vapor)"
REPORTING_PRESSURE_TRANSFORM = "P_report = observed_pressure * exp(u_pressure)"
PREDICTIVE_STATUS = "NOT_ADJUDICATED_NO_APPROVED_HELD_OUT_CUTOFF"
DIAMETER_TRANSFORM = "d_i = 3.0 angstrom + 1.0 angstrom * z_i"
BORN_RESIDUAL = "r_i = (G_i(d_i) - G_i_target) / abs(G_i_target)"
BORN_JACOBIAN = "J_ij = delta_ij * G_i_prime(d_i) * 1 angstrom / abs(G_i_target)"
WATER_FACTOR_RESIDUAL = "r_q = 1 - gamma_q_model / gamma_q_observed"
WATER_FACTOR_JACOBIAN = (
    "dr_q/df_water = -(gamma_q_model/gamma_q_observed) * dln(gamma_q_model)/df_water"
)


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
    raw_residuals: tuple[float, ...]
    scaled_residuals: tuple[float, ...]


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
    parameters: tuple[ParameterDiagnostic, ...]
    jacobian: JacobianDiagnostics
    training_rows: tuple[TrainingRowDiagnostic, ...]
    reporting_rows: tuple[ReportingRowDiagnostic, ...]
    confirmation_termination: str
    confirmation_solution_usable: bool
    confirmation_parameter_scaled_max_delta: float
    confirmation_cost_relative_delta: float
    scientific_comparison_status: str
    pressure_aad_percent: float
    liquid_density_aad_percent: float
    published_parameter_max_relative_difference: float
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
class WaterSolvationFactorRowDiagnostic:
    row_id: str
    molality_mol_per_kg: float
    observed_gamma_pm_m: float
    modeled_log_gamma_pm_m: float
    modeled_gamma_pm_m: float
    scaled_residual: float
    provider_log_derivative: float
    exact_scaled_residual_derivative: float
    reference_molality_mol_per_kg: float
    reference_convergence_error: float
    reference_derivative_convergence_error: float
    provider_fingerprint: str


@dataclass(frozen=True, slots=True)
class WaterSolvationFactorStartDiagnostic:
    name: str
    termination: str
    solution_usable: bool
    initial_cost: float
    final_cost: float
    iterations: int
    parameter: float
    rows: tuple[WaterSolvationFactorRowDiagnostic, ...]
    singular_value: float
    rank_threshold: float
    rank: int
    condition_number: float
    complete_jacobian_column: bool
    active_bound: bool
    solver_converged: bool
    failure_reasons: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class FigielWaterSolvationFactorFitResult:
    specification_id: str
    provider_fingerprint: str
    fitted_water_solvation_factor: float
    starts: tuple[WaterSolvationFactorStartDiagnostic, ...]
    start_parameter_max_abs_delta: float
    miac_rmse: float
    input_row_ids: tuple[str, ...]
    evaluated_row_ids: tuple[str, ...]
    failed_row_ids: tuple[str, ...]
    solver_converged: bool
    numerically_converged: bool
    physically_valid: bool
    workflow_valid: bool
    predictive_status: str
    failure_reasons: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class AqueousKijRowDiagnostic:
    row_id: str
    salt: str
    molality_mol_per_kg: float
    observed_gamma_pm_m: float
    modeled_log_gamma_pm_m: float
    modeled_gamma_pm_m: float
    scaled_residual: float
    local_log_derivative: tuple[float, float, float]
    reference_molality_mol_per_kg: float
    reference_convergence_error: float
    reference_derivative_convergence_error: float
    provider_fingerprint: str


@dataclass(frozen=True, slots=True)
class AqueousKijCoordinateDiagnostic:
    coordinate: int
    parameter: float
    initial_cost: float
    final_cost: float
    iterations: int
    termination: str
    solution_usable: bool
    failure_reason: str


@dataclass(frozen=True, slots=True)
class AqueousKijStartDiagnostic:
    name: str
    coordinate_order: str
    termination: str
    solution_usable: bool
    initial_cost: float
    final_cost: float
    iterations: int
    parameters: tuple[float, ...]
    coordinates: tuple[AqueousKijCoordinateDiagnostic, ...]
    rows: tuple[AqueousKijRowDiagnostic, ...]
    singular_values: tuple[float, ...]
    rank_threshold: float
    rank: int
    condition_number: float
    least_sensitive_direction: tuple[float, ...]
    complete_jacobian_columns: tuple[bool, ...]
    active_bounds: tuple[bool, ...]
    solver_converged: bool
    failure_reasons: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class AqueousKijSaltMetric:
    salt: str
    row_count: int
    miac_rmse: float
    miac_max_abs_error: float
    first_modeled_gamma_pm_m: float


@dataclass(frozen=True, slots=True)
class FigielAqueousKijFitResult:
    specification_id: str
    provider_fingerprints: tuple[str, ...]
    fitted_parameters: tuple[float, ...]
    published_parameters: tuple[float, ...]
    starts: tuple[AqueousKijStartDiagnostic, ...]
    salt_metrics: tuple[AqueousKijSaltMetric, ...]
    start_parameter_max_abs_delta: float
    published_parameter_max_abs_delta: float
    pooled_miac_rmse: float
    input_row_ids: tuple[str, ...]
    evaluated_row_ids: tuple[str, ...]
    failed_row_ids: tuple[str, ...]
    solver_converged: bool
    numerically_converged: bool
    physically_valid: bool
    workflow_valid: bool
    scientifically_valid: bool
    predictive_status: str
    recovery_status: str
    failure_reasons: tuple[str, ...]


def _born_native_payload(
    specification: BornDiameterTracerSpecification,
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
    )
    targets = tuple(
        (
            target.target_id,
            target.ion_label,
            target.active_component_id,
            target.counterion_component_id,
            target.target_j_per_mol,
            target.published_diameter_angstrom,
            target.expected_provider_fingerprint,
        )
        for target in specification.targets
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
        specification.start_diameters_angstrom,
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
    transforms = (
        (
            PARAMETER_TRANSFORM,
            "H = ((P_L-P_V)/P_obs, mu_L/RT-mu_V/RT) = 0",
            "du/dz = -H_u^{-1} H_z",
            "rho(e)=sqrt(e^2+1e-8)-1e-4 for pressure and liquid density",
            "equilibrium-newton-max-iterations=60",
            "equilibrium-line-search-max-backtracks=16",
            "equilibrium-relative-rank-threshold=1e-12",
        )
        if dataset.component_id == "monoethanolamine"
        else (
            PARAMETER_TRANSFORM,
            LIQUID_VOLUME_TRANSFORM,
            VAPOR_VOLUME_TRANSFORM,
            REPORTING_PRESSURE_TRANSFORM,
        )
    )
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
        *transforms,
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
        specification.confirmation_start,
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
        raise RuntimeError(
            "native reporting row identity did not match the immutable dataset"
        )
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
        raise RuntimeError(
            "reporting row did not belong to exactly one immutable partition"
        )
    partition = ("training", "held_out", "stress")[memberships.index(True)]
    return ReportingRowDiagnostic(
        row_id=source.row_id,
        temperature_k=source.temperature_k,
        training=source.row_id in training_ids,
        partition=partition,
        observed_pressure_pa=source.pressure_pa,
        predicted_pressure_pa=predicted_pressure,
        pressure_relative_error=(predicted_pressure - source.pressure_pa)
        / source.pressure_pa,
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
        raise ValueError(
            "model fingerprint does not match the immutable component specification"
        )
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
        raise RuntimeError(
            "compiled problem identity did not round-trip from the native solve"
        )
    variables = tuple(float(value) for value in variables_native)
    final_parameters = tuple(
        start + scale * transformed
        for start, scale, transformed in zip(
            specification.start,
            specification.parameter_scales,
            variables[: len(specification.start)],
            strict=True,
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
            liquid_molar_density_mol_m3=specification.fixed_amount_mol
            / float(native_row[3]),
            vapor_molar_density_mol_m3=specification.fixed_amount_mol
            / float(native_row[4]),
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
        if str(native_row[0]) == source.row_id
        and str(native_row[1]) == source.source_id
    )
    if training_rows_native and (
        len(training_rows_native) != len(dataset.training_rows)
        or len(training_rows) != len(dataset.training_rows)
    ):
        raise RuntimeError(
            "native training row identity did not match the immutable dataset"
        )
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
        raise RuntimeError(
            "native reporting row identity did not match the immutable dataset"
        )
    jacobian = JacobianDiagnostics(
        complete_columns=bool(complete_columns_native),
        full_singular_values=tuple(
            float(value) for value in full_singular_values_native
        ),
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
    full_columns_rank = jacobian.full_rank == len(variables)
    parameter_columns_full_rank = jacobian.parameter_rank == fitted_parameter_count
    solver_converged = (
        termination == "CONVERGENCE"
        and usable
        and math.isfinite(initial_cost)
        and math.isfinite(final_cost)
        and final_cost <= initial_cost
        and jacobian.complete_columns
        and full_columns_rank
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
    pressure_aad_percent = math.nan
    liquid_density_aad_percent = math.nan
    published_parameter_max_relative_difference = math.nan
    scientific_comparison_status = "NOT_APPLICABLE"
    if dataset.component_id == "monoethanolamine" and reporting_tuple:
        pressure_aad_percent = 100.0 * math.fsum(
            abs(row.pressure_relative_error) for row in reporting_tuple
        ) / len(reporting_tuple)
        liquid_density_aad_percent = 100.0 * math.fsum(
            abs(row.liquid_density_relative_error) for row in reporting_tuple
        ) / len(reporting_tuple)
        published = (3.0353, 3.0435, 277.174, 2586.3, 0.037470)
        published_parameter_max_relative_difference = max(
            abs(parameter.final / reference - 1.0)
            for parameter, reference in zip(parameters, published, strict=True)
        )
        scientific_comparison_status = (
            "PASSED_BOUNDED_BAYGI_INSPIRED_COMPARISON"
            if (
                published_parameter_max_relative_difference <= 0.15
                and pressure_aad_percent <= 0.62
                and liquid_density_aad_percent <= 0.12
            )
            else "FAILED_BOUNDED_BAYGI_INSPIRED_COMPARISON"
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
    if not full_columns_rank:
        failure_reasons.append(
            "training full Jacobian is rank deficient: "
            f"{jacobian.full_rank} of {len(variables)} variable columns"
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
        failure_reasons.append(
            "provider source fingerprint did not match the supplied model"
        )
        physical_valid = False
    return PureSaturationFitResult(
        component_id=dataset.component_id,
        dataset_id=dataset.dataset_id,
        specification_id=specification.specification_id,
        provider_fingerprint=observed_fingerprint,
        compiled_problem_identity=tuple(
            str(value) for value in compiled_identity_native
        ),
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
        scientific_comparison_status=scientific_comparison_status,
        pressure_aad_percent=pressure_aad_percent,
        liquid_density_aad_percent=liquid_density_aad_percent,
        published_parameter_max_relative_difference=(
            published_parameter_max_relative_difference
        ),
        failure_reasons=tuple(failure_reasons),
    )


def _born_observations(
    native_rows: tuple[object, ...],
    specification: BornDiameterTracerSpecification,
) -> tuple[BornObservationDiagnostic, ...]:
    if len(native_rows) != len(specification.targets):
        raise RuntimeError(
            "native Born observations did not match the five-target contract"
        )
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
        raise RuntimeError(
            "native Born result dimensions did not match the 5 x 5 contract"
        )
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
    if not final_cost < initial_cost:
        reasons.append("Born solve did not strictly reduce cost")
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
            raise ValueError(
                "model component order does not match the immutable Born target"
            )
        fingerprint = getattr(model, "parameter_fingerprint", None)
        if fingerprint != target.expected_provider_fingerprint:
            raise ValueError(
                "model fingerprint does not match the immutable Born target"
            )
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
    observed_fingerprints = tuple(
        row.provider_fingerprint for row in primary.observations
    )
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
        abs(row.raw_error_j_per_mol) <= specification.observable_round_trip_j_per_mol
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
        failure_reasons.append(
            "source-bound workflow identity or reference gate failed"
        )
    if not scientifically_valid:
        failure_reasons.append("source-observable reproduction gate failed")
    return BornDiameterFitResult(
        specification_id=specification.specification_id,
        compiled_problem_identity=tuple(
            str(value) for value in compiled_identity_native
        ),
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


def _figiel_aqueous_bundle(
    fixed_born_diameters_angstrom: tuple[float, ...],
    *,
    water_solvation_factor: float | None,
    bundle_id: str = "figiel-water-factor-fixed-inputs",
) -> ParameterBundle:
    catalog = ParameterBundle.from_catalog(
        "figiel-2025-reference-electrolytes", version=1
    )
    born_by_component = {
        target.active_component_id: diameter
        for target, diameter in zip(
            FIGIEL_BORN_DIAMETER_TRACER_V1.targets,
            fixed_born_diameters_angstrom,
            strict=True,
        )
    }
    records = tuple(
        replace(
            record,
            value=(born_by_component[record.component_id] * unit_registry.angstrom),
        )
        if (
            isinstance(record, SingleParameterRecord)
            and record.family == "born_diameter"
            and record.component_id in born_by_component
        )
        else (
            replace(record, value=water_solvation_factor)
            if (
                water_solvation_factor is not None
                and isinstance(record, SingleParameterRecord)
                and record.record_id == "water-solvation-factor"
            )
            else record
        )
        for record in catalog.records
    )
    return ParameterBundle.from_records(
        bundle_id=bundle_id,
        bundle_version=1,
        purpose="user-provided",
        sources=catalog.sources,
        domains=catalog.domains,
        components=catalog.components,
        singles=(
            record for record in records if isinstance(record, SingleParameterRecord)
        ),
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
                (
                    ConstantCorrelation,
                    ConstantPlusSumOfExponentialsCorrelation,
                ),
            )
        ),
        models=(
            record for record in records if isinstance(record, ModelParameterRecord)
        ),
    )


def _fixed_water_factor_model(
    specification: FigielWaterSolvationFactorSpecification,
) -> EPCSAFT:
    bundle = _figiel_aqueous_bundle(
        specification.fixed_born_diameters_angstrom,
        water_solvation_factor=None,
    )
    return EPCSAFT(bundle.select(("water", "sodium-cation", "bromide-anion")))


def _water_factor_native_payload(
    specification: FigielWaterSolvationFactorSpecification,
) -> tuple[object, ...]:
    observations = tuple(
        (
            row.row_id,
            row.molality_mol_per_kg,
            row.gamma_pm_m,
        )
        for row in specification.observations
    )
    starts = (
        ("primary", specification.starts[0]),
        ("upper", specification.starts[1]),
    )
    return (
        observations,
        specification.expected_provider_fingerprint,
        starts,
        specification.temperature_k,
        specification.pressure_pa,
        specification.parameter_bounds,
        specification.max_num_iterations,
        specification.start_wall_time_max_seconds,
        specification.function_tolerance,
        specification.gradient_tolerance,
        specification.parameter_tolerance,
        specification.rank_threshold_multiplier,
    )


def _water_factor_row(
    native_row: tuple[object, ...],
) -> WaterSolvationFactorRowDiagnostic:
    return WaterSolvationFactorRowDiagnostic(
        row_id=str(native_row[0]),
        molality_mol_per_kg=float(native_row[1]),
        observed_gamma_pm_m=float(native_row[2]),
        modeled_log_gamma_pm_m=float(native_row[3]),
        modeled_gamma_pm_m=float(native_row[4]),
        scaled_residual=float(native_row[5]),
        provider_log_derivative=float(native_row[6]),
        exact_scaled_residual_derivative=float(native_row[7]),
        reference_molality_mol_per_kg=float(native_row[8]),
        reference_convergence_error=float(native_row[9]),
        reference_derivative_convergence_error=float(native_row[10]),
        provider_fingerprint=str(native_row[11]),
    )


def _water_factor_start(
    native_start: tuple[object, ...],
    specification: FigielWaterSolvationFactorSpecification,
) -> WaterSolvationFactorStartDiagnostic:
    if len(native_start) != 13:
        raise RuntimeError("native water-factor start has the wrong dimension")
    rows = tuple(_water_factor_row(tuple(row)) for row in tuple(native_start[7]))
    parameter = float(native_start[6])
    initial_cost = float(native_start[3])
    final_cost = float(native_start[4])
    singular_value = float(native_start[8])
    rank_threshold = float(native_start[9])
    termination = str(native_start[1])
    solution_usable = bool(native_start[2])
    rank = int(native_start[10])
    complete_column = bool(native_start[11])
    active_tolerance = math.sqrt(math.ulp(1.0)) * 2.0
    active_bound = (
        min(
            parameter - specification.parameter_bounds[0],
            specification.parameter_bounds[1] - parameter,
        )
        <= active_tolerance
    )
    finite = all(
        math.isfinite(value)
        for value in (
            parameter,
            initial_cost,
            final_cost,
            singular_value,
            rank_threshold,
            *(
                value
                for row in rows
                for value in (
                    row.modeled_log_gamma_pm_m,
                    row.modeled_gamma_pm_m,
                    row.scaled_residual,
                    row.provider_log_derivative,
                    row.exact_scaled_residual_derivative,
                    row.reference_molality_mol_per_kg,
                    row.reference_convergence_error,
                    row.reference_derivative_convergence_error,
                )
            ),
        )
    )
    failure_reasons: list[str] = []
    if termination != "CONVERGENCE":
        failure_reasons.append(f"Ceres termination was {termination}")
    if not solution_usable:
        failure_reasons.append("Ceres solution was unusable")
    if not finite:
        failure_reasons.append("water-factor solution or Jacobian was nonfinite")
    if final_cost > initial_cost + math.ulp(max(1.0, abs(initial_cost))):
        failure_reasons.append("water-factor solve increased cost")
    if len(rows) != 21:
        failure_reasons.append("water-factor result did not contain 21 rows")
    if not complete_column:
        failure_reasons.append("water-factor Jacobian column was incomplete")
    if rank != 1:
        failure_reasons.append(f"water-factor Jacobian rank was {rank} of 1")
    native_failure = str(native_start[12]).strip()
    if native_failure:
        failure_reasons.append(native_failure)
    return WaterSolvationFactorStartDiagnostic(
        name=str(native_start[0]),
        termination=termination,
        solution_usable=solution_usable,
        initial_cost=initial_cost,
        final_cost=final_cost,
        iterations=int(native_start[5]),
        parameter=parameter,
        rows=rows,
        singular_value=singular_value,
        rank_threshold=rank_threshold,
        rank=rank,
        condition_number=1.0 if rank == 1 else math.inf,
        complete_jacobian_column=complete_column,
        active_bound=active_bound,
        solver_converged=(termination == "CONVERGENCE" and solution_usable),
        failure_reasons=tuple(failure_reasons),
    )


def fit_figiel_water_solvation_factor() -> FigielWaterSolvationFactorFitResult:
    specification = FIGIEL_WATER_SOLVATION_FACTOR_V1
    model = _fixed_water_factor_model(specification)
    expected_fingerprint = specification.expected_provider_fingerprint
    if model.parameter_fingerprint != expected_fingerprint:
        raise RuntimeError("installed Provider model fingerprint does not match")
    payload = _water_factor_native_payload(specification)
    native_starts = _native.solve_figiel_water_factor(native_sdk(model), payload)
    starts = tuple(
        _water_factor_start(tuple(start), specification)
        for start in tuple(native_starts)
    )
    if tuple(start.name for start in starts) != ("primary", "upper"):
        raise RuntimeError("native water-factor starts did not match the contract")

    primary = starts[0]
    start_delta = abs(starts[1].parameter - primary.parameter)
    input_row_ids = tuple(row.row_id for row in specification.observations)
    evaluated_row_ids = tuple(row.row_id for row in primary.rows if row.row_id)
    failed_row_ids = tuple(
        row_id for row_id in input_row_ids if row_id not in evaluated_row_ids
    )
    solver_converged = all(start.solver_converged for start in starts)
    numerically_converged = (
        solver_converged
        and all(not start.failure_reasons for start in starts)
        and start_delta <= specification.start_agreement_max_abs
    )
    physically_valid = all(
        row.modeled_gamma_pm_m > 0.0
        and row.provider_fingerprint == expected_fingerprint
        and math.isfinite(row.reference_molality_mol_per_kg)
        and math.isfinite(row.reference_convergence_error)
        and math.isfinite(row.reference_derivative_convergence_error)
        for start in starts
        for row in start.rows
    )
    workflow_valid = (
        input_row_ids == evaluated_row_ids
        and not failed_row_ids
        and len(primary.rows) == 21
    )
    miac_rmse = math.sqrt(
        sum(
            (row.modeled_gamma_pm_m - row.observed_gamma_pm_m) ** 2
            for row in primary.rows
        )
        / len(primary.rows)
    )
    failure_reasons = [
        f"{start.name}: {reason}"
        for start in starts
        for reason in start.failure_reasons
    ]
    if start_delta > specification.start_agreement_max_abs:
        failure_reasons.append("declared-start agreement gate failed")
    if not physically_valid:
        failure_reasons.append("Provider state or reference diagnostic gate failed")
    if not workflow_valid:
        failure_reasons.append("source-bound identity or row accounting failed")
    return FigielWaterSolvationFactorFitResult(
        specification_id=specification.specification_id,
        provider_fingerprint=expected_fingerprint,
        fitted_water_solvation_factor=primary.parameter,
        starts=starts,
        start_parameter_max_abs_delta=start_delta,
        miac_rmse=miac_rmse,
        input_row_ids=input_row_ids,
        evaluated_row_ids=evaluated_row_ids,
        failed_row_ids=failed_row_ids,
        solver_converged=solver_converged,
        numerically_converged=numerically_converged,
        physically_valid=physically_valid,
        workflow_valid=workflow_valid,
        predictive_status=PREDICTIVE_STATUS,
        failure_reasons=tuple(failure_reasons),
    )


def _aqueous_kij_models(
    specification: FigielAqueousKijSpecification,
) -> tuple[EPCSAFT, ...]:
    bundle = _figiel_aqueous_bundle(
        specification.fixed_born_diameters_angstrom,
        water_solvation_factor=specification.fixed_water_solvation_factor,
        bundle_id="figiel-stage-c-fixed-inputs",
    )
    return tuple(
        EPCSAFT(bundle.select(("water", cation, anion)))
        for _, cation, anion, _ in specification.salt_contracts
    )


def _aqueous_kij_native_payload(
    specification: FigielAqueousKijSpecification,
) -> tuple[object, ...]:
    salt_indices = {
        salt: index
        for index, (salt, _, _, _) in enumerate(specification.salt_contracts)
    }
    columns = {
        salt: coordinate_indices
        for salt, _, _, coordinate_indices in specification.salt_contracts
    }
    observations = tuple(
        (
            row.row_id,
            row.salt,
            salt_indices[row.salt],
            row.molality_mol_per_kg,
            row.gamma_pm_m,
            columns[row.salt],
        )
        for row in specification.observations
    )
    schedules = tuple(
        (
            name,
            (value,) * len(specification.coordinate_order),
            order,
        )
        for name, value, order in specification.start_schedules
    )
    return (
        observations,
        specification.expected_provider_fingerprints,
        schedules,
        specification.temperature_k,
        specification.pressure_pa,
        specification.parameter_bounds,
        specification.published_parameters,
        specification.max_num_iterations,
        specification.start_wall_time_max_seconds,
        specification.function_tolerance,
        specification.gradient_tolerance,
        specification.parameter_tolerance,
        specification.rank_threshold_multiplier,
        specification.specification_id,
    )


def _aqueous_kij_row(
    native_row: tuple[object, ...],
) -> AqueousKijRowDiagnostic:
    local_derivative = tuple(float(value) for value in native_row[7])
    if len(local_derivative) != 3:
        raise RuntimeError("native aqueous-kij derivative has the wrong dimension")
    return AqueousKijRowDiagnostic(
        row_id=str(native_row[0]),
        salt=str(native_row[1]),
        molality_mol_per_kg=float(native_row[2]),
        observed_gamma_pm_m=float(native_row[3]),
        modeled_log_gamma_pm_m=float(native_row[4]),
        modeled_gamma_pm_m=float(native_row[5]),
        scaled_residual=float(native_row[6]),
        local_log_derivative=(
            local_derivative[0],
            local_derivative[1],
            local_derivative[2],
        ),
        reference_molality_mol_per_kg=float(native_row[8]),
        reference_convergence_error=float(native_row[9]),
        reference_derivative_convergence_error=float(native_row[10]),
        provider_fingerprint=str(native_row[11]),
    )


def _aqueous_kij_start(
    native_start: tuple[object, ...],
    specification: FigielAqueousKijSpecification,
) -> AqueousKijStartDiagnostic:
    if len(native_start) != 17:
        raise RuntimeError("native aqueous-kij start has the wrong dimension")
    parameters = tuple(float(value) for value in native_start[7])
    rows = tuple(_aqueous_kij_row(tuple(row)) for row in tuple(native_start[8]))
    singular_values = tuple(float(value) for value in native_start[9])
    least_sensitive = tuple(float(value) for value in native_start[13])
    complete_columns = tuple(bool(value) for value in native_start[14])
    coordinates = tuple(
        AqueousKijCoordinateDiagnostic(
            coordinate=int(item[0]),
            parameter=float(item[1]),
            initial_cost=float(item[2]),
            final_cost=float(item[3]),
            iterations=int(item[4]),
            termination=str(item[5]),
            solution_usable=bool(item[6]),
            failure_reason=str(item[7]),
        )
        for item in (tuple(value) for value in tuple(native_start[15]))
    )
    if not (
        len(parameters)
        == len(singular_values)
        == len(least_sensitive)
        == len(complete_columns)
        == 11
        and len(coordinates) == 11
        and tuple(item.coordinate for item in coordinates) == tuple(range(11))
        and len(rows) == 164
    ):
        raise RuntimeError(
            "native aqueous-kij result dimensions did not match 164 x 11"
        )
    termination = str(native_start[2])
    solution_usable = bool(native_start[3])
    initial_cost = float(native_start[4])
    final_cost = float(native_start[5])
    rank_threshold = float(native_start[10])
    rank = int(native_start[11])
    condition_number = float(native_start[12])
    bound_tolerance = math.sqrt(math.ulp(1.0)) * 2.0
    active_bounds = tuple(
        min(
            parameter - specification.parameter_bounds[0],
            specification.parameter_bounds[1] - parameter,
        )
        <= bound_tolerance
        for parameter in parameters
    )
    finite = all(
        math.isfinite(value)
        for value in (
            initial_cost,
            final_cost,
            rank_threshold,
            condition_number,
            *parameters,
            *singular_values,
            *least_sensitive,
            *(
                value
                for coordinate in coordinates
                for value in (
                    coordinate.parameter,
                    coordinate.initial_cost,
                    coordinate.final_cost,
                )
            ),
            *(
                value
                for row in rows
                for value in (
                    row.modeled_log_gamma_pm_m,
                    row.modeled_gamma_pm_m,
                    row.scaled_residual,
                    *row.local_log_derivative,
                    row.reference_molality_mol_per_kg,
                    row.reference_convergence_error,
                    row.reference_derivative_convergence_error,
                )
            ),
        )
    )
    failure_reasons: list[str] = []
    if termination != "CONVERGENCE":
        failure_reasons.append(f"Ceres termination was {termination}")
    if not solution_usable:
        failure_reasons.append("Ceres solution was unusable")
    if not finite:
        failure_reasons.append("aqueous-kij solution or Jacobian was nonfinite")
    if final_cost > initial_cost + math.ulp(max(1.0, abs(initial_cost))):
        failure_reasons.append("aqueous-kij solve increased cost")
    for coordinate in coordinates:
        if coordinate.termination != "CONVERGENCE":
            failure_reasons.append(
                f"coordinate {coordinate.coordinate} terminated "
                f"{coordinate.termination}"
            )
        if not coordinate.solution_usable:
            failure_reasons.append(
                f"coordinate {coordinate.coordinate} solution was unusable"
            )
        if coordinate.failure_reason:
            failure_reasons.append(
                f"coordinate {coordinate.coordinate}: {coordinate.failure_reason}"
            )
        if coordinate.final_cost > coordinate.initial_cost + math.ulp(
            max(1.0, abs(coordinate.initial_cost))
        ):
            failure_reasons.append(f"coordinate {coordinate.coordinate} increased cost")
        if not (
            specification.parameter_bounds[0]
            <= coordinate.parameter
            <= specification.parameter_bounds[1]
        ):
            failure_reasons.append(
                f"coordinate {coordinate.coordinate} left its bounds"
            )
    if not all(complete_columns):
        failure_reasons.append("aqueous-kij Jacobian columns were incomplete")
    if rank != 11:
        failure_reasons.append(f"aqueous-kij Jacobian rank was {rank} of 11")
    native_failure = str(native_start[16]).strip()
    if native_failure:
        failure_reasons.append(native_failure)
    return AqueousKijStartDiagnostic(
        name=str(native_start[0]),
        coordinate_order=str(native_start[1]),
        termination=termination,
        solution_usable=solution_usable,
        initial_cost=initial_cost,
        final_cost=final_cost,
        iterations=int(native_start[6]),
        parameters=parameters,
        coordinates=coordinates,
        rows=rows,
        singular_values=singular_values,
        rank_threshold=rank_threshold,
        rank=rank,
        condition_number=condition_number,
        least_sensitive_direction=least_sensitive,
        complete_jacobian_columns=complete_columns,
        active_bounds=active_bounds,
        solver_converged=(termination == "CONVERGENCE" and solution_usable),
        failure_reasons=tuple(failure_reasons),
    )


def fit_figiel_aqueous_kij() -> FigielAqueousKijFitResult:
    specification = FIGIEL_AQUEOUS_KIJ_V1
    models = _aqueous_kij_models(specification)
    fingerprints = tuple(model.parameter_fingerprint for model in models)
    if fingerprints != specification.expected_provider_fingerprints:
        raise RuntimeError("installed Provider model fingerprints do not match Stage C")
    payload = _aqueous_kij_native_payload(specification)
    capsules = tuple(native_sdk(model) for model in models)
    primary = _aqueous_kij_start(
        tuple(_native.solve_figiel_kij(capsules, payload, 0)),
        specification,
    )
    confirmation = _aqueous_kij_start(
        tuple(_native.solve_figiel_kij(capsules, payload, 1)),
        specification,
    )
    starts = (primary, confirmation)
    if tuple(start.name for start in starts) != ("primary", "confirmation"):
        raise RuntimeError("native aqueous-kij starts did not match the contract")

    input_row_ids = tuple(row.row_id for row in specification.observations)
    evaluated_row_ids = tuple(row.row_id for row in primary.rows)
    failed_row_ids = tuple(
        row_id for row_id in input_row_ids if row_id not in evaluated_row_ids
    )
    start_delta = max(
        abs(value - reference)
        for start in starts[1:]
        for value, reference in zip(start.parameters, primary.parameters, strict=True)
    )
    published_delta = max(
        abs(value - published)
        for value, published in zip(
            primary.parameters,
            specification.published_parameters,
            strict=True,
        )
    )
    rows_by_salt = {
        salt: tuple(row for row in primary.rows if row.salt == salt)
        for salt, *_ in specification.salt_contracts
    }
    complete_reporting_rows = (
        input_row_ids == evaluated_row_ids
        and not failed_row_ids
        and all(rows_by_salt.values())
    )
    salt_metrics = tuple(
        AqueousKijSaltMetric(
            salt=salt,
            row_count=len(salt_rows),
            miac_rmse=math.sqrt(
                sum(
                    (row.modeled_gamma_pm_m - row.observed_gamma_pm_m) ** 2
                    for row in salt_rows
                )
                / len(salt_rows)
            ),
            miac_max_abs_error=max(
                abs(row.modeled_gamma_pm_m - row.observed_gamma_pm_m)
                for row in salt_rows
            ),
            first_modeled_gamma_pm_m=salt_rows[0].modeled_gamma_pm_m,
        )
        for salt, *_ in specification.salt_contracts
        for salt_rows in (rows_by_salt[salt],)
        if salt_rows
    )
    pooled_rmse = (
        math.sqrt(
            sum(
                (row.modeled_gamma_pm_m - row.observed_gamma_pm_m) ** 2
                for row in primary.rows
            )
            / len(primary.rows)
        )
        if complete_reporting_rows
        else math.inf
    )
    solver_converged = all(start.solver_converged for start in starts)
    numerically_converged = (
        solver_converged
        and all(not start.failure_reasons for start in starts)
        and start_delta <= specification.start_agreement_max_abs
    )
    physically_valid = all(
        row.modeled_gamma_pm_m > 0.0
        and row.provider_fingerprint
        == fingerprints[
            next(
                index
                for index, contract in enumerate(specification.salt_contracts)
                if contract[0] == row.salt
            )
        ]
        for start in starts
        for row in start.rows
    )
    workflow_valid = complete_reporting_rows and len(primary.rows) == 164
    observable_gates = (
        pooled_rmse <= 0.17
        and all(metric.miac_rmse <= 0.35 for metric in salt_metrics)
        and all(metric.miac_max_abs_error <= 1.25 for metric in salt_metrics)
        and all(metric.first_modeled_gamma_pm_m < 0.98 for metric in salt_metrics)
    )
    scientific_gates = (
        published_delta <= specification.published_parameter_max_abs_delta
        and observable_gates
    )
    scientifically_valid = (
        numerically_converged
        and physically_valid
        and workflow_valid
        and scientific_gates
    )
    failure_reasons = [
        f"{start.name}: {reason}"
        for start in starts
        for reason in start.failure_reasons
    ]
    if start_delta > specification.start_agreement_max_abs:
        failure_reasons.append("declared-start agreement gate failed")
    if not physically_valid:
        failure_reasons.append("Provider state or fingerprint gate failed")
    if not workflow_valid:
        failure_reasons.append("source-bound row accounting failed")
    if published_delta > specification.published_parameter_max_abs_delta:
        failure_reasons.append(
            "fitted interactions did not reproduce the printed Table 4/5 tuple"
        )
    if not observable_gates:
        failure_reasons.append("source-observable reproduction gate failed")
    return FigielAqueousKijFitResult(
        specification_id=specification.specification_id,
        provider_fingerprints=fingerprints,
        fitted_parameters=primary.parameters,
        published_parameters=specification.published_parameters,
        starts=starts,
        salt_metrics=salt_metrics,
        start_parameter_max_abs_delta=start_delta,
        published_parameter_max_abs_delta=published_delta,
        pooled_miac_rmse=pooled_rmse,
        input_row_ids=input_row_ids,
        evaluated_row_ids=evaluated_row_ids,
        failed_row_ids=failed_row_ids,
        solver_converged=solver_converged,
        numerically_converged=numerically_converged,
        physically_valid=physically_valid,
        workflow_valid=workflow_valid,
        scientifically_valid=scientifically_valid,
        predictive_status=PREDICTIVE_STATUS,
        recovery_status=(
            "FIGIEL_AQUEOUS_KIJ_CONDITIONALLY_RECOVERED"
            if scientifically_valid
            else ("SOURCE_DESCRIBED_STAGED_RECOVERY_DID_NOT_REPRODUCE_PRINTED_TUPLE")
        ),
        failure_reasons=tuple(failure_reasons),
    )
