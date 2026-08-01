from __future__ import annotations

import ctypes
import hashlib
import json
from collections import Counter
from pathlib import Path

import pytest
from epcsaft import Mixture, Parameters, native_sdk
from parameter_cases import aqueous_parameters

import epcsaft_regression
from epcsaft_regression import _native
from epcsaft_regression.workflow import _aqueous_kij_native_payload


def _models() -> tuple[Mixture, ...]:
    return tuple(
        Mixture(
            Parameters.from_dictionary(
                aqueous_parameters(("water", cation, anion))
            )
        )
        for _, cation, anion, _ in epcsaft_regression.FIGIEL_AQUEOUS_KIJ_V1.salt_contracts
    )


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
        ("evaluation_budget_size", ctypes.c_size_t),
        ("evaluate_aqueous_miac_kij_batch_bounded", ctypes.c_void_p),
        (
            "evaluate_aqueous_miac_solvation_factor_batch_bounded",
            ctypes.c_void_p,
        ),
    ]


def _copied_provider_capsules(
    capsules: tuple[object, ...],
    *,
    clear_scalar: bool = False,
    truncate_before_bounded_batch: bool = False,
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
        assert table.table_size >= ctypes.sizeof(table)
        table.table_size = ctypes.sizeof(table)
        if clear_scalar:
            table.evaluate_aqueous_miac_kij = None
        if truncate_before_bounded_batch:
            table.table_size = _NativeSdkV1.evaluate_aqueous_miac_kij_batch_bounded.offset
        tables.append(table)
        copied_capsules.append(new_capsule(ctypes.addressof(table), name, None))
    return tuple(copied_capsules), tuple(tables)


def test_aqueous_kij_contract_is_eleven_parameters_over_all_rows() -> None:
    specification = epcsaft_regression.FIGIEL_AQUEOUS_KIJ_V1

    assert specification.specification_id == "figiel-2025-aqueous-kij-v1"
    assert specification.source_validation_commit == (
        "8944d34f7002cda1bb8760e606cc1f11696f58cd"
    )
    assert specification.source_hamer_wu_csv_sha256 == (
        "2f63e13f06a5b0f4e8bca2980b6a8d9d7fb0f839153c43e3a71952daf9796595"
    )
    assert specification.fixed_born_diameters_angstrom == (
        2.7888130173797934,
        3.4524616464076425,
        4.147266741279482,
        4.101505615791675,
        4.476998527506598,
    )
    assert specification.fixed_water_solvation_factor == 1.5590515389548207
    assert specification.fixed_born_evidence_subject_sha256 == (
        "55ea2cd69af62c45b26179cfab6939760de23058b5a7e8c880a79f67faa417ed"
    )
    assert specification.fixed_water_factor_regression_commit == (
        "882e0735ed1b5586a591682da1fd3d78f46636d4"
    )
    assert len(specification.observations) == 164
    assert Counter(row.salt for row in specification.observations) == {
        "LiCl": 29,
        "NaCl": 29,
        "KCl": 28,
        "LiBr": 29,
        "NaBr": 21,
        "KBr": 28,
    }
    assert specification.coordinate_order == (
        ("water", "lithium-cation"),
        ("water", "sodium-cation"),
        ("water", "potassium-cation"),
        ("water", "chloride-anion"),
        ("water", "bromide-anion"),
        ("lithium-cation", "chloride-anion"),
        ("sodium-cation", "chloride-anion"),
        ("potassium-cation", "chloride-anion"),
        ("lithium-cation", "bromide-anion"),
        ("sodium-cation", "bromide-anion"),
        ("potassium-cation", "bromide-anion"),
    )
    assert specification.published_parameters == (
        -0.4,
        -0.3,
        -0.1,
        -0.3,
        -0.3,
        0.8,
        0.8,
        0.0,
        0.5,
        0.65,
        -0.35,
    )
    assert specification.parameter_bounds == (-1.0, 1.0)
    assert specification.start_schedules == (
        ("primary", 0.0, "forward"),
        ("confirmation", 0.25, "reverse"),
    )
    assert specification.start_agreement_max_abs == 1.0e-5
    assert specification.published_parameter_max_abs_delta == 0.05
    assert specification.expected_provider_fingerprints == (
        "sha256:00d5bc709841bffd8cbde9df6afafaa9520d23a500be59aae54ecd63b4369ad9",
        "sha256:569f4ace0d776e8b631df0ff563a2686356b28ad5e24556c38909b371fd811ac",
        "sha256:17955d601749bc94802325dbc88cdc62b3ee954f1e5e1ad02c0547d8187b27d4",
        "sha256:31f40e1441358b8729c0963c808cb5290cde6d455f6f910a5c7abe95315a8aea",
        "sha256:8fedd3d395a20591b05ee7fff2c1e698872c97b07e82a0668923d656f818790e",
        "sha256:c0af15d984b8324ebc4b622a9ab68153c43eb13026bd74a6bcd08bba7d49dea5",
    )


def test_aqueous_kij_runtime_is_one_private_native_owner() -> None:
    assert hasattr(epcsaft_regression, "fit_figiel_aqueous_kij")
    assert not hasattr(epcsaft_regression, "persist_provider_parameters")
    assert not hasattr(epcsaft_regression, "fit_generic_parameters")
    fingerprints = ("sha256:" + "0" * 64,) * 6
    assert _aqueous_kij_native_payload(
        epcsaft_regression.FIGIEL_AQUEOUS_KIJ_V1, fingerprints
    )[-1] == "figiel-2025-aqueous-kij-v1"


@pytest.mark.campaign
def test_aqueous_kij_requires_only_the_bounded_batch_callback() -> None:
    specification = epcsaft_regression.FIGIEL_AQUEOUS_KIJ_V1
    models = _models()
    capsules = tuple(native_sdk(model) for model in models)
    payload = _aqueous_kij_native_payload(
        specification, tuple(model.parameter_fingerprint for model in models)
    )
    copied_capsules, keepalive = _copied_provider_capsules(
        capsules, clear_scalar=True
    )

    residuals, jacobian, rows = _native.evaluate_figiel_kij(
        copied_capsules, payload, specification.published_parameters
    )

    assert keepalive
    assert len(residuals) == 164
    assert len(jacobian) == 164 * 11
    assert len(rows) == 164


def test_aqueous_kij_rejects_provider_truncated_before_bounded_batch() -> None:
    specification = epcsaft_regression.FIGIEL_AQUEOUS_KIJ_V1
    models = _models()
    capsules = tuple(native_sdk(model) for model in models)
    payload = _aqueous_kij_native_payload(
        specification, tuple(model.parameter_fingerprint for model in models)
    )
    copied_capsules, keepalive = _copied_provider_capsules(
        capsules, truncate_before_bounded_batch=True
    )

    with pytest.raises(
        RuntimeError, match="provider aqueous-kij capability unavailable"
    ):
        _native.evaluate_figiel_kij(
            copied_capsules, payload, specification.published_parameters
        )
    assert keepalive


@pytest.mark.campaign
def test_aqueous_kij_residual_jacobian_maps_all_exact_provider_columns() -> None:
    specification = epcsaft_regression.FIGIEL_AQUEOUS_KIJ_V1
    models = _models()
    capsules = tuple(native_sdk(model) for model in models)
    payload = _aqueous_kij_native_payload(
        specification, tuple(model.parameter_fingerprint for model in models)
    )

    _, jacobian_native, rows_native = _native.evaluate_figiel_kij(
        capsules, payload, specification.published_parameters
    )
    jacobian = tuple(float(value) for value in jacobian_native)
    column_map = {salt: columns for salt, _, _, columns in specification.salt_contracts}
    for row_index, row_native in enumerate(rows_native):
        modeled = float(row_native[5])
        observed = float(row_native[3])
        local = tuple(float(value) for value in row_native[7])
        expected = [0.0] * 11
        for column, derivative in zip(
            column_map[str(row_native[1])], local, strict=True
        ):
            expected[column] = -(modeled / observed) * derivative
        actual = jacobian[row_index * 11 : (row_index + 1) * 11]
        assert actual == pytest.approx(expected, rel=1.0e-14, abs=1.0e-14)


@pytest.mark.campaign
def test_aqueous_kij_exact_jacobian_matches_callback_value_direction() -> None:
    specification = epcsaft_regression.FIGIEL_AQUEOUS_KIJ_V1
    models = _models()
    capsules = tuple(native_sdk(model) for model in models)
    payload = _aqueous_kij_native_payload(
        specification, tuple(model.parameter_fingerprint for model in models)
    )
    parameters = specification.published_parameters
    direction = tuple(
        (-1.0 if index % 2 else 1.0) * (index + 1.0) / 11.0
        for index in range(11)
    )
    _, jacobian_native, _ = _native.evaluate_figiel_kij(
        capsules, payload, parameters
    )
    jacobian = tuple(float(value) for value in jacobian_native)
    differences: list[tuple[float, ...]] = []
    for step in (1.0e-4, 5.0e-5):
        plus = tuple(
            value + step * component
            for value, component in zip(parameters, direction, strict=True)
        )
        minus = tuple(
            value - step * component
            for value, component in zip(parameters, direction, strict=True)
        )
        residuals_plus, _, _ = _native.evaluate_figiel_kij(capsules, payload, plus)
        residuals_minus, _, _ = _native.evaluate_figiel_kij(
            capsules, payload, minus
        )
        differences.append(
            tuple(
                (float(plus_value) - float(minus_value)) / (2.0 * step)
                for plus_value, minus_value in zip(
                    residuals_plus, residuals_minus, strict=True
                )
            )
        )
    for row in range(164):
        exact = sum(
            jacobian[row * 11 + column] * direction[column]
            for column in range(11)
        )
        coarse = differences[0][row]
        fine = differences[1][row]
        tolerance = max(
            1.0e-8,
            20.0 * abs(coarse - fine),
            2.0e-8 * abs(exact),
        )
        assert abs(exact - fine) <= tolerance


@pytest.mark.campaign
def test_public_aqueous_kij_fit_replays_retained_negative_result() -> None:
    result = epcsaft_regression.fit_figiel_aqueous_kij(models=_models())

    assert result.solver_converged
    assert result.numerically_converged
    assert result.physically_valid
    assert result.workflow_valid
    assert not result.scientifically_valid
    assert result.recovery_status == (
        "SOURCE_DESCRIBED_STAGED_RECOVERY_DID_NOT_REPRODUCE_PRINTED_TUPLE"
    )
    assert result.input_row_ids == result.evaluated_row_ids
    assert result.failed_row_ids == ()
    assert tuple(start.rank for start in result.starts) == (11, 11)
    assert result.start_parameter_max_abs_delta <= 1.0e-5
    assert result.published_parameter_max_abs_delta > 0.05


def test_conditional_campaign_evidence_retains_separate_statuses() -> None:
    evidence_path = (
        Path(__file__).parents[1]
        / "evidence"
        / "figiel-aqueous-kij-conditional-recovery.json"
    )
    evidence_bytes = evidence_path.read_bytes()
    assert hashlib.sha256(evidence_bytes).hexdigest() == (
        "448a3d48e661a2a40e01833147a3b0842edd4436663b953cf2c90dde178d0417"
    )
    evidence = json.loads(evidence_bytes)

    assert evidence["primary"]["rank"] == 11
    assert evidence["primary"]["complete_jacobian_columns"] == [True] * 11
    assert evidence["primary"]["active_bounds"] == [
        False,
        False,
        False,
        False,
        False,
        True,
        False,
        False,
        False,
        False,
        False,
    ]
    assert evidence["problem"]["fixed_born_diameters_angstrom"] == list(
        epcsaft_regression.FIGIEL_AQUEOUS_KIJ_V1.fixed_born_diameters_angstrom
    )
    assert evidence["problem"]["fixed_water_solvation_factor"] == (
        epcsaft_regression.FIGIEL_AQUEOUS_KIJ_V1.fixed_water_solvation_factor
    )
    assert evidence["derivative_checks"]["passed_entries"] == 1804
    assert evidence["derivative_checks"]["total_entries"] == 1804
    assert evidence["regression_artifact"]["implementation_commit"] == (
        "846c3671a437f961f92c20300c4cb3c874ff6009"
    )
    assert evidence["regression_artifact"]["wheel_sha256"] == (
        "0fe3ddfe6a31993c7859609d2bdc13cb6b69bd037b276671cbee2b1ded32a365"
    )
    assert evidence["comparison"]["start_max_abs_delta"] <= 1.0e-5
    assert evidence["comparison"]["published_max_abs_delta"] > 0.05
    assert evidence["comparison"]["coordinates_within_gate"] == 2
    assert evidence["observables"]["evaluated_rows"] == 164
    assert evidence["statuses"] == {
        "solver_converged": True,
        "numerically_converged": True,
        "physically_valid": True,
        "workflow_valid": True,
        "scientifically_valid": False,
        "predictive_status": "NOT_ADJUDICATED_NO_APPROVED_HELD_OUT_CUTOFF",
    }
