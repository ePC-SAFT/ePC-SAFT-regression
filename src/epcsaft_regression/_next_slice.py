from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path


class NextSliceReadiness(StrEnum):
    READY_WAITING_PROPANE_PRIMARY_SOURCE_PACKET = (
        "READY_WAITING_PROPANE_PRIMARY_SOURCE_PACKET"
    )
    READY_WAITING_PROVIDER_CORRECTION = "READY_WAITING_PROVIDER_CORRECTION"


@dataclass(frozen=True, slots=True)
class NextSlicePreflight:
    propane_status: NextSliceReadiness
    binary_kij_status: NextSliceReadiness
    provider_coordinate_order: tuple[str, ...]


def _callback_declaration(header: str, callback_type: str) -> str:
    callback_start = header.find(f"typedef int (*{callback_type})(")
    if callback_start < 0:
        return ""
    callback_end = header.find(");", callback_start)
    if callback_end < 0:
        return ""
    return header[callback_start:callback_end]


def _contains_parameters_in_order(
    declaration: str, parameters: tuple[str, ...]
) -> bool:
    offset = 0
    for parameter in parameters:
        offset = declaration.find(parameter, offset)
        if offset < 0:
            return False
        offset += len(parameter)
    return True


def _provider_coordinate_order(header: str) -> tuple[str, ...]:
    active_kij_callback = _callback_declaration(
        header, "epcsaft_evaluate_mixture_phase_kij_v1"
    )
    if _contains_parameters_in_order(
        active_kij_callback,
        (
            "const double* amounts_mol",
            "size_t amount_count",
            "double volume_m3",
            "double k_ij",
            "epcsaft_mixture_phase_block_result_v1* result",
        ),
    ):
        return (
            "methane_amount_mol",
            "ethane_amount_mol",
            "volume_m3",
            "kij",
        )

    fixed_callback = _callback_declaration(
        header, "epcsaft_evaluate_mixture_phase_v1"
    )
    if _contains_parameters_in_order(
        fixed_callback,
        (
            "const double* amounts_mol",
            "size_t amount_count",
            "double volume_m3",
            "epcsaft_mixture_phase_block_result_v1* result",
        ),
    ):
        return (
            "component_amounts_mol_in_parameter_order",
            "volume_m3",
        )
    return ()


def preflight_next_slice(*, provider_header: Path) -> NextSlicePreflight:
    coordinate_order = _provider_coordinate_order(
        provider_header.read_text(encoding="utf-8")
    )
    return NextSlicePreflight(
        propane_status=(
            NextSliceReadiness.READY_WAITING_PROPANE_PRIMARY_SOURCE_PACKET
        ),
        binary_kij_status=NextSliceReadiness.READY_WAITING_PROVIDER_CORRECTION,
        provider_coordinate_order=coordinate_order,
    )
