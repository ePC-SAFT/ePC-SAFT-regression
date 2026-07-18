from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
import hashlib
from pathlib import Path


class NextSliceReadiness(StrEnum):
    READY_WAITING_PROPANE_PRIMARY_SOURCE_PACKET = (
        "READY_WAITING_PROPANE_PRIMARY_SOURCE_PACKET"
    )
    READY_WAITING_PROVIDER_ACTIVE_KIJ_CALLBACK = (
        "READY_WAITING_PROVIDER_ACTIVE_KIJ_CALLBACK"
    )
    READY_FOR_PROPANE_IMPLEMENTATION = "READY_FOR_PROPANE_IMPLEMENTATION"
    READY_FOR_BINARY_KIJ_IMPLEMENTATION = "READY_FOR_BINARY_KIJ_IMPLEMENTATION"


@dataclass(frozen=True, slots=True)
class BinaryKijFitSpecification:
    specification_id: str
    dataset_id: str
    source_owner_commit: str
    source_csv_sha256: str
    source_packet_sha256: str
    source_pair_value: float
    source_pair_record_sha256: str
    component_order: tuple[str, str]
    expected_provider_fingerprint: str
    provider_coordinate_order: tuple[str, str, str, str]
    row_count: int
    variables_per_row: int
    residuals_per_row: int
    derivative_order: int
    kij_start: float
    kij_bounds: tuple[float, float]
    kij_scale: float
    liquid_volume_reference_m3: float
    vapor_volume_reference: str
    liquid_volume_bounds_m3: tuple[float, float]
    vapor_volume_bounds_m3: tuple[float, float]
    residual_names: tuple[str, str, str, str]
    residual_scales: tuple[str, str, str, str]
    residual_weights: tuple[float, float, float, float]
    required_full_rank: int
    required_projected_parameter_rank: int
    rank_threshold: str
    kij_bound_margin: float
    topology_relative_separation_min: float
    pressure_scaled_residual_max: float
    confirmation_volume_start_multipliers: tuple[float, float]
    confirmation_kij_starts: tuple[float, float]
    confirmation_kij_scaled_max_delta: float
    confirmation_cost_relative_delta: float
    max_num_iterations: int
    function_tolerance: float
    gradient_tolerance: float
    parameter_tolerance: float
    ceres_linear_solver: str
    ceres_num_threads: int
    predictive_status: str

    @property
    def variable_count(self) -> int:
        return 1 + self.variables_per_row * self.row_count

    @property
    def residual_count(self) -> int:
        return self.residuals_per_row * self.row_count


BINARY_KIJ_FIT_V1 = BinaryKijFitSpecification(
    specification_id="methane-ethane-constant-kij-lifted-volumes-v1",
    dataset_id="may-2015-methane-ethane-vle",
    source_owner_commit="73a37f5935e919a34d1e4fa3af285951d6fac8e7",
    source_csv_sha256="5cd1e74925a3c6504f5106dcf911f2cae2d6e99a5133fccc20454d8991bdbc7f",
    source_packet_sha256="d43433e93b354e01f96d330c760818a24b775026461ce795e45774cfb11ac94e",
    source_pair_value=0.0,
    source_pair_record_sha256="747e8281c7a1e4240ee4badbc0bedd047521fb303726699b38c99fccf7f74c2a",
    component_order=("methane", "ethane"),
    expected_provider_fingerprint=(
        "sha256:307fcb28d535b94782f3e3caf4012c0c8c0dc87ee4239d6c316de56553543286"
    ),
    provider_coordinate_order=(
        "methane_amount_mol",
        "ethane_amount_mol",
        "volume_m3",
        "kij",
    ),
    row_count=17,
    variables_per_row=2,
    residuals_per_row=4,
    derivative_order=2,
    kij_start=0.0,
    kij_bounds=(-0.15, 0.10),
    kij_scale=0.01,
    liquid_volume_reference_m3=6.5e-5,
    vapor_volume_reference="R*T_K/P_Pa",
    liquid_volume_bounds_m3=(2.0e-5, 1.0e-4),
    vapor_volume_bounds_m3=(1.0e-4, 1.0e-2),
    residual_names=(
        "liquid_pressure_closure",
        "vapor_pressure_closure",
        "methane_mu_over_rt_equality",
        "ethane_mu_over_rt_equality",
    ),
    residual_scales=("observed_pressure_pa", "observed_pressure_pa", "1", "1"),
    residual_weights=(0.25, 0.25, 0.25, 0.25),
    required_full_rank=35,
    required_projected_parameter_rank=1,
    rank_threshold="sigma_max*max(residual_count,variable_count)*epsilon_double",
    kij_bound_margin=1.0e-8,
    topology_relative_separation_min=1.0e-3,
    pressure_scaled_residual_max=1.0e-8,
    confirmation_volume_start_multipliers=(1.01, 0.98),
    confirmation_kij_starts=(-0.05, 0.05),
    confirmation_kij_scaled_max_delta=1.0e-5,
    confirmation_cost_relative_delta=1.0e-8,
    max_num_iterations=500,
    function_tolerance=1.0e-10,
    gradient_tolerance=1.0e-10,
    parameter_tolerance=1.0e-10,
    ceres_linear_solver="DENSE_QR",
    ceres_num_threads=1,
    predictive_status="NOT_ADJUDICATED_NO_APPROVED_HELD_OUT_CUTOFF",
)


@dataclass(frozen=True, slots=True)
class BinaryKijJacobianDiagnostics:
    complete_columns: bool
    full_singular_values: tuple[float, ...]
    full_rank: int
    full_condition_number: float
    projected_parameter_singular_value: float
    projected_parameter_rank: int


@dataclass(frozen=True, slots=True)
class BinaryKijRowDiagnostic:
    row_id: str
    temperature_k: float
    observed_pressure_pa: float
    liquid_methane_mole_fraction: float
    vapor_methane_mole_fraction: float
    liquid_volume_m3: float
    vapor_volume_m3: float
    liquid_pressure_pa: float
    vapor_pressure_pa: float
    liquid_mu_over_rt: tuple[float, float]
    vapor_mu_over_rt: tuple[float, float]
    liquid_stability_slope: float
    vapor_stability_slope: float
    raw_residuals: tuple[float, float, float, float]
    scaled_residuals: tuple[float, float, float, float]
    physically_valid: bool
    failure_reasons: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class BinaryKijFitResult:
    dataset_id: str
    specification_id: str
    source_csv_sha256: str
    source_packet_sha256: str
    provider_fingerprint: str
    provider_header_sha256: str
    provider_wheel_sha256: str
    solver_converged: bool
    numerically_converged: bool
    physically_valid: bool
    predictive_status: str
    termination: str
    solution_usable: bool
    initial_cost: float
    final_cost: float
    iterations: int
    kij_start: float
    kij_final: float
    kij_movement: float
    kij_bounds: tuple[float, float]
    kij_active_bound: str | None
    jacobian: BinaryKijJacobianDiagnostics
    rows: tuple[BinaryKijRowDiagnostic, ...]
    confirmation_termination: tuple[str, ...]
    confirmation_solution_usable: bool
    confirmation_kij_scaled_max_delta: float
    confirmation_cost_relative_delta: float
    failure_reasons: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class NextSliceUpstreamEvidence:
    propane_validation_commit: str = ""
    propane_primary_source_packet_sha256: str = ""
    propane_source_csv_sha256: str = ""
    provider_active_kij_commit: str = ""
    provider_active_kij_header_sha256: str = ""
    provider_active_kij_wheel_sha256: str = ""


@dataclass(frozen=True, slots=True)
class NextSlicePreflight:
    propane_status: NextSliceReadiness
    binary_kij_status: NextSliceReadiness
    provider_header_sha256: str
    provider_coordinate_order: tuple[str, ...]
    missing_upstream_hashes: tuple[str, ...]


def _provider_coordinate_order(header: str) -> tuple[str, ...]:
    callback_start = header.find("typedef int (*epcsaft_evaluate_mixture_phase_v1)(")
    if callback_start < 0:
        return ()
    callback_end = header.find(");", callback_start)
    if callback_end < 0:
        return ()
    callback = header[callback_start:callback_end]
    if "double k_ij" in callback:
        return BINARY_KIJ_FIT_V1.provider_coordinate_order
    return ("component_amounts_mol_in_parameter_order", "volume_m3")


def preflight_next_slice(
    *,
    provider_header: Path,
    evidence: NextSliceUpstreamEvidence = NextSliceUpstreamEvidence(),
) -> NextSlicePreflight:
    header_bytes = provider_header.read_bytes()
    coordinate_order = _provider_coordinate_order(header_bytes.decode("utf-8"))
    supplied = {
        field: getattr(evidence, field)
        for field in NextSliceUpstreamEvidence.__dataclass_fields__
    }
    missing = tuple(field for field, value in supplied.items() if not value)
    propane_ready = not any(field.startswith("propane_") for field in missing)
    binary_ready = (
        coordinate_order == BINARY_KIJ_FIT_V1.provider_coordinate_order
        and not any(field.startswith("provider_active_kij_") for field in missing)
    )
    return NextSlicePreflight(
        propane_status=(
            NextSliceReadiness.READY_FOR_PROPANE_IMPLEMENTATION
            if propane_ready
            else NextSliceReadiness.READY_WAITING_PROPANE_PRIMARY_SOURCE_PACKET
        ),
        binary_kij_status=(
            NextSliceReadiness.READY_FOR_BINARY_KIJ_IMPLEMENTATION
            if binary_ready
            else NextSliceReadiness.READY_WAITING_PROVIDER_ACTIVE_KIJ_CALLBACK
        ),
        provider_header_sha256=hashlib.sha256(header_bytes).hexdigest(),
        provider_coordinate_order=coordinate_order,
        missing_upstream_hashes=missing,
    )
