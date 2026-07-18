from __future__ import annotations

from dataclasses import FrozenInstanceError
from pathlib import Path

import epcsaft
import pytest

from epcsaft_regression._next_slice import (
    BINARY_KIJ_FIT_V1,
    BinaryKijFitResult,
    BinaryKijJacobianDiagnostics,
    BinaryKijRowDiagnostic,
    NextSliceReadiness,
    preflight_next_slice,
)


def test_binary_kij_contract_freezes_the_approved_problem() -> None:
    specification = BINARY_KIJ_FIT_V1

    assert specification.specification_id == "methane-ethane-constant-kij-lifted-volumes-v1"
    assert specification.dataset_id == "may-2015-methane-ethane-vle"
    assert specification.source_owner_commit == "73a37f5935e919a34d1e4fa3af285951d6fac8e7"
    assert specification.source_csv_sha256 == (
        "5cd1e74925a3c6504f5106dcf911f2cae2d6e99a5133fccc20454d8991bdbc7f"
    )
    assert specification.source_packet_sha256 == (
        "d43433e93b354e01f96d330c760818a24b775026461ce795e45774cfb11ac94e"
    )
    assert specification.component_order == ("methane", "ethane")
    assert specification.expected_provider_fingerprint == (
        "sha256:307fcb28d535b94782f3e3caf4012c0c8c0dc87ee4239d6c316de56553543286"
    )
    assert specification.provider_coordinate_order == (
        "methane_amount_mol",
        "ethane_amount_mol",
        "volume_m3",
        "kij",
    )
    assert specification.row_count == 17
    assert specification.variables_per_row == 2
    assert specification.residuals_per_row == 4
    assert specification.variable_count == 35
    assert specification.residual_count == 68
    assert specification.derivative_order == 2
    assert specification.kij_start == 0.0
    assert specification.kij_bounds == (-0.15, 0.10)
    assert specification.kij_scale == 0.01
    assert specification.liquid_volume_reference_m3 == 6.5e-5
    assert specification.vapor_volume_reference == "R*T_K/P_Pa"
    assert specification.residual_names == (
        "liquid_pressure_closure",
        "vapor_pressure_closure",
        "methane_mu_over_rt_equality",
        "ethane_mu_over_rt_equality",
    )
    assert specification.residual_scales == (
        "observed_pressure_pa",
        "observed_pressure_pa",
        "1",
        "1",
    )
    assert specification.residual_weights == (0.25, 0.25, 0.25, 0.25)
    assert specification.required_full_rank == 35
    assert specification.required_projected_parameter_rank == 1
    assert specification.confirmation_volume_start_multipliers == (1.01, 0.98)
    assert specification.confirmation_kij_starts == (-0.05, 0.05)
    assert specification.confirmation_kij_scaled_max_delta == 1.0e-5
    assert specification.confirmation_cost_relative_delta == 1.0e-8

    with pytest.raises(FrozenInstanceError):
        specification.kij_start = 0.01  # type: ignore[misc]

    for contract_type in (
        BinaryKijJacobianDiagnostics,
        BinaryKijRowDiagnostic,
        BinaryKijFitResult,
    ):
        assert contract_type.__dataclass_params__.frozen


def test_current_upstreams_report_exact_waiting_gates() -> None:
    provider_header = (
        Path(epcsaft.__file__).parent / "include" / "epcsaft" / "native_sdk_v1.h"
    )

    result = preflight_next_slice(provider_header=provider_header)

    assert result.propane_status is NextSliceReadiness.READY_WAITING_PROPANE_PRIMARY_SOURCE_PACKET
    assert result.binary_kij_status is NextSliceReadiness.READY_WAITING_PROVIDER_ACTIVE_KIJ_CALLBACK
    assert result.provider_header_sha256 == (
        "b1dc4a666799603ea87fe515ff08226455a11e52908b7109ebe1841369cb92df"
    )
    assert result.provider_coordinate_order == (
        "component_amounts_mol_in_parameter_order",
        "volume_m3",
    )
    assert result.missing_upstream_hashes == (
        "propane_validation_commit",
        "propane_primary_source_packet_sha256",
        "propane_source_csv_sha256",
        "provider_active_kij_commit",
        "provider_active_kij_header_sha256",
        "provider_active_kij_wheel_sha256",
    )
