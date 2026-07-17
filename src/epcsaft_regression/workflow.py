from __future__ import annotations

from dataclasses import dataclass
import math

from epcsaft import native_sdk

from . import _native
from .records import (
    EXPECTED_PACKAGED_DATA_SHA256,
    MethaneFitSpecification,
    MethaneSaturationDataset,
)


PROVIDER_COMMIT = "4b10cb899c94687cae734980285badb224dc95e6"
PROVIDER_WHEEL_SHA256 = "f92f79c8d6f614660e5c201b7061c9b02b5cd1a25a4ed8c8fee0b59adaabf2bf"
PROVIDER_CAPSULE = "epcsaft.native_sdk.v1"
PARAMETER_TRANSFORM = "p_j = start_j + parameter_scale_j * z_j"
LIQUID_VOLUME_TRANSFORM = "V_liquid = (molar_mass / observed_liquid_density) * exp(u_liquid)"
VAPOR_VOLUME_TRANSFORM = "V_vapor = (R * T / observed_pressure) * exp(u_vapor)"
REPORTING_PRESSURE_TRANSFORM = "P_report = observed_pressure * exp(u_pressure)"


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
        PROVIDER_COMMIT,
        PROVIDER_WHEEL_SHA256,
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
        tuple(
            (
                row.row_id,
                row.species,
                row.temperature_k,
                row.pressure_pa,
                row.liquid_density_kg_m3,
                row.source_id,
            )
            for row in dataset.training_rows
        ),
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
    reporting_payload = tuple(
        (
            row.row_id,
            row.species,
            row.temperature_k,
            row.pressure_pa,
            row.liquid_density_kg_m3,
            row.source_id,
        )
        for row in dataset.rows
    )
    transport = _native.solve(capsule, payload, reporting_payload)
    if tuple(transport[23]) != payload[0]:
        raise RuntimeError("compiled problem identity did not round-trip from the native solve")
    variables = tuple(float(value) for value in transport[5])
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
        for source, native_row in zip(dataset.training_rows, transport[8], strict=True)
        if str(native_row[0]) == source.row_id and str(native_row[1]) == source.source_id
    )
    if len(training_rows) != len(dataset.training_rows):
        raise RuntimeError("native training row identity did not match the immutable dataset")
    training_ids = frozenset(dataset.training_row_ids)
    reporting_rows: list[ReportingRowDiagnostic] = []
    for source, native_row in zip(dataset.rows, transport[20], strict=True):
        reasons: list[str] = []
        if str(native_row[0]) != source.row_id or str(native_row[1]) != source.source_id:
            raise RuntimeError("native reporting row identity did not match the immutable dataset")
        termination = str(native_row[12])
        usable = bool(native_row[13])
        liquid_volume = float(native_row[7])
        vapor_volume = float(native_row[8])
        liquid_slope = float(native_row[9])
        vapor_slope = float(native_row[10])
        raw = tuple(float(value) for value in native_row[11])
        if termination != "CONVERGENCE":
            reasons.append(f"reporting Ceres termination was {termination}")
        if not usable:
            reasons.append("reporting Ceres solution was unusable")
        if not liquid_volume < vapor_volume:
            reasons.append("reporting phases lost volume ordering")
        if liquid_slope <= 0.0 or vapor_slope <= 0.0:
            reasons.append("reporting phase was mechanically unstable")
        if (
            abs(raw[0]) / source.pressure_pa
            > specification.reporting_pressure_scaled_residual_max
            or abs(raw[1]) / source.pressure_pa
            > specification.reporting_pressure_scaled_residual_max
        ):
            reasons.append("reporting scaled pressure closure exceeded its threshold")
        if abs(raw[2]) > specification.reporting_chemical_potential_residual_max:
            reasons.append("reporting chemical-potential closure exceeded its threshold")
        if not all(math.isfinite(value) for value in (*raw, *native_row[2:11])):
            reasons.append("reporting diagnostics were nonfinite")
        predicted_pressure = float(native_row[5])
        predicted_density = float(native_row[6])
        reporting_rows.append(
            ReportingRowDiagnostic(
                row_id=source.row_id,
                temperature_k=source.temperature_k,
                training=source.row_id in training_ids,
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
                liquid_molar_density_mol_m3=specification.fixed_amount_mol / liquid_volume,
                vapor_molar_density_mol_m3=specification.fixed_amount_mol / vapor_volume,
                liquid_mass_density_kg_m3=specification.methane_molar_mass_kg_per_mol
                / liquid_volume,
                vapor_mass_density_kg_m3=specification.methane_molar_mass_kg_per_mol
                / vapor_volume,
                liquid_stability_slope=liquid_slope,
                vapor_stability_slope=vapor_slope,
                raw_equilibrium_residuals=raw,
                termination=termination,
                solution_usable=usable,
                physically_valid=not reasons,
                failure_reasons=tuple(reasons),
            )
        )
    jacobian = JacobianDiagnostics(
        complete_columns=bool(transport[15]),
        full_singular_values=tuple(float(value) for value in transport[9]),
        full_rank=int(transport[10]),
        full_condition_number=float(transport[11]),
        parameter_singular_values=tuple(float(value) for value in transport[12]),
        parameter_rank=int(transport[13]),
        parameter_condition_number=float(transport[14]),
    )
    termination = str(transport[0])
    usable = bool(transport[1])
    initial_cost = float(transport[2])
    final_cost = float(transport[3])
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
    )
    confirmation_termination = str(transport[18])
    confirmation_usable = bool(transport[19])
    parameter_delta = float(transport[16])
    cost_delta = float(transport[17])
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
    if not solver_converged:
        failure_reasons.append("training solver convergence gate failed")
    if not numerically_converged:
        failure_reasons.append("confirmation solve numerical convergence gate failed")
    if not physical_valid:
        failure_reasons.append("training or reporting physical validity gate failed")
    if str(transport[21]) != getattr(model, "parameter_fingerprint", None):
        failure_reasons.append("provider source fingerprint did not match the supplied model")
        physical_valid = False
    return MethaneFitResult(
        dataset_id=dataset.dataset_id,
        specification_id=specification.specification_id,
        provider_fingerprint=str(transport[21]),
        compiled_problem_identity=tuple(str(value) for value in transport[23]),
        solver_converged=solver_converged,
        numerically_converged=numerically_converged,
        physically_valid=physical_valid,
        termination=termination,
        solution_usable=usable,
        initial_cost=initial_cost,
        final_cost=final_cost,
        iterations=int(transport[4]),
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
