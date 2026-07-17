from __future__ import annotations

from dataclasses import dataclass
import math

from epcsaft import native_sdk

from . import _native
from .records import (
    EXPECTED_PACKAGED_DATA_SHA256,
    MethaneFitSpecification,
    MethaneSaturationDataset,
    SaturationObservation,
)


PROVIDER_CAPSULE = "epcsaft.native_sdk.v1"
PARAMETER_TRANSFORM = "p_j = start_j + parameter_scale_j * z_j"
LIQUID_VOLUME_TRANSFORM = "V_liquid = (molar_mass / observed_liquid_density) * exp(u_liquid)"
VAPOR_VOLUME_TRANSFORM = "V_vapor = (R * T / observed_pressure) * exp(u_vapor)"
REPORTING_PRESSURE_TRANSFORM = "P_report = observed_pressure * exp(u_pressure)"


def _row_payload(row: SaturationObservation) -> tuple[object, ...]:
    return (
        row.row_id,
        row.species,
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
class MethaneFitResult:
    dataset_id: str
    specification_id: str
    provider_fingerprint: str
    compiled_problem_identity: tuple[str, ...]
    solver_converged: bool
    numerically_converged: bool
    physically_valid: bool
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


def _native_payload(
    dataset: MethaneSaturationDataset,
    specification: MethaneFitSpecification,
    provider_fingerprint: str,
) -> tuple[object, ...]:
    identity = (
        dataset.dataset_id,
        dataset.species,
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
        EXPECTED_PACKAGED_DATA_SHA256,
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
        specification.methane_molar_mass_kg_per_mol,
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
    specification: MethaneFitSpecification,
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
    return ReportingRowDiagnostic(
        row_id=source.row_id,
        temperature_k=source.temperature_k,
        training=source.row_id in training_ids,
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
            specification.methane_molar_mass_kg_per_mol * liquid_molar_density
        ),
        vapor_mass_density_kg_m3=(
            specification.methane_molar_mass_kg_per_mol * vapor_molar_density
        ),
        liquid_stability_slope=liquid_slope,
        vapor_stability_slope=vapor_slope,
        raw_equilibrium_residuals=raw,
        termination=termination,
        solution_usable=usable,
        physically_valid=not reasons,
        failure_reasons=tuple(reasons),
    )


def fit_methane_saturation(
    *,
    model: object,
    dataset: MethaneSaturationDataset,
    specification: MethaneFitSpecification,
) -> MethaneFitResult:
    if type(dataset) is not MethaneSaturationDataset:
        raise TypeError("dataset must be an exact MethaneSaturationDataset")
    if type(specification) is not MethaneFitSpecification:
        raise TypeError("specification must be an exact MethaneFitSpecification")
    capsule = native_sdk(model)
    provider_fingerprint = getattr(model, "parameter_fingerprint", None)
    if not isinstance(provider_fingerprint, str) or not provider_fingerprint:
        raise ValueError("model must expose a nonblank provider parameter_fingerprint")
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
            liquid_mass_density_kg_m3=specification.methane_molar_mass_kg_per_mol
            / float(native_row[3]),
            vapor_mass_density_kg_m3=specification.methane_molar_mass_kg_per_mol
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
    reporting_rows = [
        _reporting_row_diagnostic(source, training_ids, specification, native_row)
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
    solver_converged = (
        termination == "CONVERGENCE"
        and usable
        and math.isfinite(initial_cost)
        and math.isfinite(final_cost)
        and final_cost <= initial_cost
        and jacobian.complete_columns
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
    physical_valid = (
        solver_converged
        and all(row.liquid_volume_m3 < row.vapor_volume_m3 for row in training_rows)
        and all(row.liquid_stability_slope > 0.0 for row in training_rows)
        and all(row.vapor_stability_slope > 0.0 for row in training_rows)
        and all(row.physically_valid for row in reporting_tuple)
    )
    failure_reasons: list[str] = []
    if native_failure_reason:
        failure_reasons.append(native_failure_reason)
    failure_reasons.extend(
        f"{row.row_id}: {reason}"
        for row in reporting_tuple
        for reason in row.failure_reasons
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
    return MethaneFitResult(
        dataset_id=dataset.dataset_id,
        specification_id=specification.specification_id,
        provider_fingerprint=observed_fingerprint,
        compiled_problem_identity=tuple(str(value) for value in compiled_identity_native),
        solver_converged=solver_converged,
        numerically_converged=numerically_converged,
        physically_valid=physical_valid,
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
