from __future__ import annotations

from collections import Counter
import ctypes
from dataclasses import replace
import hashlib
from importlib.resources import files
from pathlib import Path
import sys

import pytest

import epcsaft_regression
import epcsaft_regression._native as native
from epcsaft import native_sdk
from epcsaft_regression.records import (
    FIGIEL_AQUEOUS_KIJ_COORDINATES,
    FIGIEL_AQUEOUS_MIAC_PACKAGED_SHA256,
    FIGIEL_AQUEOUS_PUBLISHED_KIJ,
    FIGIEL_STAGED_AQUEOUS_RECOVERY_V1,
)
from epcsaft_regression.workflow import (
    AQUEOUS_KIJ_COLUMNS,
    _aqueous_native_payload,
    _figiel_models,
)

TOOLS_DIRECTORY = Path(__file__).parents[1] / "tools"
sys.path.insert(0, str(TOOLS_DIRECTORY))
from run_figiel_born_candidate import _staged_conclusion  # noqa: E402
sys.path.pop(0)


class _NativeSdkV1(ctypes.Structure):
    _fields_ = [
        ("abi_version", ctypes.c_uint32),
        ("table_size", ctypes.c_size_t),
        ("result_size", ctypes.c_size_t),
        ("model_context", ctypes.c_void_p),
        ("evaluate_pure_phase", ctypes.c_void_p),
        ("parameterized_result_size", ctypes.c_size_t),
        ("evaluate_pure_phase_parameters", ctypes.c_void_p),
        ("component_count", ctypes.c_size_t),
        ("mixture_result_size", ctypes.c_size_t),
        ("evaluate_mixture_phase", ctypes.c_void_p),
        ("evaluate_mixture_phase_kij", ctypes.c_void_p),
        ("component_ids", ctypes.c_void_p),
        ("component_charges", ctypes.c_void_p),
        ("evaluate_electrolyte_phase", ctypes.c_void_p),
        ("evaluate_molar_volume_bounds", ctypes.c_void_p),
        ("evaluate_packing_fraction", ctypes.c_void_p),
        ("source_temperature_min_k", ctypes.c_double),
        ("source_temperature_max_k", ctypes.c_double),
        ("total_ion_mole_fraction_max", ctypes.c_double),
        ("ion_solvation_born_result_size", ctypes.c_size_t),
        ("evaluate_ion_solvation_born", ctypes.c_void_p),
        ("aqueous_miac_kij_result_size", ctypes.c_size_t),
        ("evaluate_aqueous_miac_kij", ctypes.c_void_p),
        ("neutral_reference_basis_row_count", ctypes.c_size_t),
        ("neutral_reference_result_size", ctypes.c_size_t),
        ("evaluate_neutral_reference", ctypes.c_void_p),
        ("aqueous_miac_solvation_factor_result_size", ctypes.c_size_t),
        ("evaluate_aqueous_miac_solvation_factor", ctypes.c_void_p),
        ("evaluate_aqueous_miac_kij_batch", ctypes.c_void_p),
        ("evaluate_aqueous_miac_solvation_factor_batch", ctypes.c_void_p),
    ]


def _copied_provider_capsules(
    capsules: tuple[object, ...],
    *,
    clear_scalar: str | None = None,
    truncate_before: str | None = None,
) -> tuple[tuple[object, ...], tuple[_NativeSdkV1, ...]]:
    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    new_capsule = ctypes.pythonapi.PyCapsule_New
    new_capsule.argtypes = (ctypes.c_void_p, ctypes.c_char_p, ctypes.c_void_p)
    new_capsule.restype = ctypes.py_object
    name = b"epcsaft.native_sdk.v1"
    tables: list[_NativeSdkV1] = []
    copied_capsules: list[object] = []
    for capsule in capsules:
        pointer = get_pointer(capsule, name)
        table = _NativeSdkV1()
        ctypes.memmove(ctypes.addressof(table), pointer, ctypes.sizeof(table))
        assert table.table_size == ctypes.sizeof(table)
        if clear_scalar is not None:
            setattr(table, clear_scalar, None)
        if truncate_before is not None:
            table.table_size = getattr(_NativeSdkV1, truncate_before).offset
        tables.append(table)
        copied_capsules.append(
            new_capsule(ctypes.addressof(table), name, None)
        )
    return tuple(copied_capsules), tuple(tables)


def test_staged_figiel_contract_freezes_all_rows_coordinates_and_controls() -> None:
    specification = FIGIEL_STAGED_AQUEOUS_RECOVERY_V1
    packaged = files("epcsaft_regression.data").joinpath(
        "figiel-aqueous-miac-targets.csv"
    ).read_bytes()

    assert hashlib.sha256(packaged).hexdigest() == FIGIEL_AQUEOUS_MIAC_PACKAGED_SHA256
    assert specification.source_ledger_sha256 == (
        "f405a3e48d21cd979a8dd480d5f8cb3be40754f5d6babf368b505b5f305607f0"
    )
    assert specification.source_parameter_packet_sha256 == (
        "932e8baa90fcefbaa8c3a8730cdeadd83a4c01f0a3b109f4e4cd0319aee9312b"
    )
    assert specification.source_metadata_sha256 == (
        "8ea06c6ca5452d01448a03f9a76cf7d0c35bb99c9abe23ccb1729d56c71d468f"
    )
    assert specification.source_si_extraction_sha256 == (
        "85bd39f727158d5a9d6eea6828c1673f73850e783a655b09660cc9b66d84321a"
    )
    assert len(specification.observations) == 164
    assert len({row.row_id for row in specification.observations}) == 164
    assert Counter(row.salt for row in specification.observations) == {
        "LiCl": 29,
        "NaCl": 29,
        "KCl": 28,
        "LiBr": 29,
        "NaBr": 21,
        "KBr": 28,
    }
    assert len(specification.stage_b_observations) == 21
    assert specification.kij_coordinate_order == FIGIEL_AQUEOUS_KIJ_COORDINATES
    assert specification.published_kij == FIGIEL_AQUEOUS_PUBLISHED_KIJ
    assert specification.solvent_factor_bounds == (1.0, 2.0)
    assert specification.solvent_factor_starts == (1.2, 1.8)
    assert specification.kij_bounds == (-1.0, 1.0)
    assert specification.kij_starts == (
        (0.0,) * 11,
        (-0.5,) * 11,
        (0.5,) * 11,
    )
    assert specification.max_confirmation_cycles == 3
    assert specification.cycle_scaled_max_delta == 1.0e-5
    assert specification.aqueous_start_wall_time_max_seconds == 180.0
    assert all(row.gamma_pm_m > 0.0 for row in specification.observations)


def test_staged_figiel_contract_rejects_row_or_gate_mutation() -> None:
    specification = FIGIEL_STAGED_AQUEOUS_RECOVERY_V1
    with pytest.raises(ValueError, match="164 unique rows"):
        replace(specification, observations=specification.observations[:-1])
    with pytest.raises(ValueError, match="frozen design"):
        replace(specification, parameter_comparison_max_abs=0.05001)
    with pytest.raises(ValueError, match="frozen design"):
        replace(specification, aqueous_start_wall_time_max_seconds=181.0)


def test_public_surface_is_one_closed_staged_workflow() -> None:
    assert (
        epcsaft_regression.FIGIEL_STAGED_AQUEOUS_RECOVERY_V1
        is FIGIEL_STAGED_AQUEOUS_RECOVERY_V1
    )
    assert hasattr(epcsaft_regression, "FigielStagedAqueousRecoveryResult")
    assert hasattr(epcsaft_regression, "fit_figiel_staged_aqueous_parameters")
    for excluded in (
        "ParameterRegistry",
        "ParameterOverlay",
        "fit_electrolyte_parameters",
        "persist_provider_parameters",
    ):
        assert not hasattr(epcsaft_regression, excluded)


def test_staged_evidence_conclusion_fails_closed_on_formula_checks() -> None:
    class Result:
        solver_converged = True
        numerically_converged = True
        physically_valid = True
        workflow_valid = True
        scientifically_valid = True

    assert _staged_conclusion(Result(), [{"passed": True}]) == (
        "FIGIEL_STAGED_AQUEOUS_RECOVERY_CANDIDATE_PASSED"
    )
    assert _staged_conclusion(Result(), [{"passed": False}]) == (
        "BLOCKED_STRICT_LOCAL_GATES"
    )


def test_staged_transient_models_consume_exact_provider_aqueous_jacobians() -> None:
    specification = FIGIEL_STAGED_AQUEOUS_RECOVERY_V1
    born = tuple(
        target.published_diameter_angstrom
        for target in epcsaft_regression.FIGIEL_BORN_DIAMETER_TRACER_V1.targets
    )
    _, models = _figiel_models(born, 1.5, specification.published_kij)
    fingerprints = tuple(model.parameter_fingerprint for model in models)
    provider_capsules = tuple(native_sdk(model) for model in models)

    for stage, parameters, starts, scalar_field in (
        (
            "solvation_factor",
            (1.5,),
            (("primary", (1.5,)),),
            "evaluate_aqueous_miac_solvation_factor",
        ),
        (
            "aqueous_kij",
            specification.published_kij,
            (("primary", (0.0,) * 11),),
            "evaluate_aqueous_miac_kij",
        ),
    ):
        capsules, tables = _copied_provider_capsules(
            provider_capsules, clear_scalar=scalar_field
        )
        payload = _aqueous_native_payload(
            specification,
            stage=stage,
            expected_fingerprints=fingerprints,
            starts=starts,
        )
        residuals, jacobian, rows, identity = native.evaluate_figiel_aqueous(
            capsules, payload, parameters
        )
        parameter_count = 1 if stage == "solvation_factor" else 11
        expected_rows = 21 if stage == "solvation_factor" else 164
        assert tuple(identity) == payload[0]
        assert len(residuals) == len(rows) == expected_rows
        assert len(jacobian) == expected_rows * parameter_count
        for row_index, row in enumerate(rows):
            ratio = float(row[5]) / float(row[3])
            assert float(residuals[row_index]) == 1.0 - ratio
            columns = (0,) if stage == "solvation_factor" else AQUEOUS_KIJ_COLUMNS[str(row[1])]
            expected = [0.0] * parameter_count
            for column, derivative in zip(columns, row[7], strict=True):
                expected[column] = -ratio * float(derivative)
            observed_jacobian_row = tuple(
                jacobian[
                    row_index * parameter_count:(row_index + 1) * parameter_count
                ]
            )
            assert observed_jacobian_row == tuple(expected)
        assert tables


@pytest.mark.parametrize(
    ("stage", "parameters", "truncate_before"),
    (
        (
            "solvation_factor",
            (1.5,),
            "evaluate_aqueous_miac_solvation_factor_batch",
        ),
        (
            "aqueous_kij",
            FIGIEL_AQUEOUS_PUBLISHED_KIJ,
            "evaluate_aqueous_miac_kij_batch",
        ),
    ),
)
def test_staged_aqueous_fit_rejects_provider_without_required_batch_tail(
    stage: str,
    parameters: tuple[float, ...],
    truncate_before: str,
) -> None:
    specification = FIGIEL_STAGED_AQUEOUS_RECOVERY_V1
    born = tuple(
        target.published_diameter_angstrom
        for target in epcsaft_regression.FIGIEL_BORN_DIAMETER_TRACER_V1.targets
    )
    _, models = _figiel_models(born, 1.5, specification.published_kij)
    fingerprints = tuple(model.parameter_fingerprint for model in models)
    provider_capsules = tuple(native_sdk(model) for model in models)
    capsules, tables = _copied_provider_capsules(
        provider_capsules, truncate_before=truncate_before
    )
    payload = _aqueous_native_payload(
        specification,
        stage=stage,
        expected_fingerprints=fingerprints,
        starts=(("primary", parameters),),
    )

    with pytest.raises(RuntimeError, match="provider capability unavailable"):
        native.evaluate_figiel_aqueous(capsules, payload, parameters)
    assert tables
