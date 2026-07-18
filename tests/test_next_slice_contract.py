from __future__ import annotations

from pathlib import Path

import epcsaft
import pytest

from epcsaft_regression._next_slice import (
    NextSliceReadiness,
    preflight_next_slice,
)


PROVIDER_REPOSITORY_HEADER = (
    Path(__file__).resolve().parents[2]
    / "ePC-SAFT"
    / "src"
    / "epcsaft"
    / "include"
    / "epcsaft"
    / "native_sdk_v1.h"
)


def test_installed_provider_remains_correction_gated() -> None:
    provider_header = (
        Path(epcsaft.__file__).parent / "include" / "epcsaft" / "native_sdk_v1.h"
    )

    result = preflight_next_slice(provider_header=provider_header)

    assert result.propane_status is NextSliceReadiness.READY_WAITING_PROPANE_PRIMARY_SOURCE_PACKET
    assert result.binary_kij_status is NextSliceReadiness.READY_WAITING_PROVIDER_CORRECTION
    assert result.provider_coordinate_order == (
        "component_amounts_mol_in_parameter_order",
        "volume_m3",
    )


@pytest.mark.skipif(
    not PROVIDER_REPOSITORY_HEADER.is_file(),
    reason="the bounded cross-repository provider header is unavailable",
)
def test_real_provider_header_detects_the_separate_active_kij_callback() -> None:
    result = preflight_next_slice(provider_header=PROVIDER_REPOSITORY_HEADER)

    assert result.binary_kij_status is NextSliceReadiness.READY_WAITING_PROVIDER_CORRECTION
    assert result.provider_coordinate_order == (
        "methane_amount_mol",
        "ethane_amount_mol",
        "volume_m3",
        "kij",
    )
