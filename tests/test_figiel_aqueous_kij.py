from __future__ import annotations

from collections import Counter
import hashlib
import json
from pathlib import Path

from epcsaft import native_sdk
import epcsaft_regression
import pytest
from epcsaft_regression import _native
from epcsaft_regression.workflow import (
    _aqueous_kij_models,
    _aqueous_kij_native_payload,
)


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
        "sha256:2e86bfea7ba8d860482189c50b5c8c9ab20736e216d7dff511fc32fc8b66156f",
        "sha256:7a88365c8a238c0d3d650a6d9b49e477e9f5e762e51abf779d0f65fe3aa1cb73",
        "sha256:f598133ee5894859c95c58ff03c0fc9c8be48b6033890d26a753c2db3c94d2fe",
        "sha256:0f3795c2e83813a7c06d2d12d534284e729b5b2c836bfe9602ab6f3be4806807",
        "sha256:1102deed6b08dbab82db24aea4187056db4432b9674187fb126cf19470ed65e1",
        "sha256:bef74bbf5556b32a1e440da9500c4e9740af457d370ec5ee900c0e59bcd79a25",
    )


def test_aqueous_kij_runtime_is_one_private_native_owner() -> None:
    assert hasattr(epcsaft_regression, "fit_figiel_aqueous_kij")
    assert not hasattr(epcsaft_regression, "persist_provider_parameters")
    assert not hasattr(epcsaft_regression, "fit_generic_parameters")


def test_libr_interaction_fit_is_start_confirmed_but_misses_printed_value() -> None:
    specification = epcsaft_regression.FIGIEL_AQUEOUS_KIJ_V1
    models = _aqueous_kij_models(specification)
    capsules = tuple(native_sdk(model) for model in models)
    payload = _aqueous_kij_native_payload(specification)

    primary = _native.solve_figiel_kij_coordinate(capsules, payload, 0, 8)
    confirmation = _native.solve_figiel_kij_coordinate(capsules, payload, 1, 8)

    assert primary[5:8] == ("CONVERGENCE", True, "")
    assert confirmation[5:8] == ("CONVERGENCE", True, "")
    assert primary[1] == pytest.approx(0.7607943631938054, rel=0.0, abs=1.0e-9)
    assert confirmation[1] == pytest.approx(0.7607929822737367, rel=0.0, abs=1.0e-9)
    assert abs(primary[1] - confirmation[1]) <= 1.0e-5
    assert abs(primary[1] - specification.published_parameters[8]) > 0.05


def test_aqueous_kij_residual_jacobian_maps_all_exact_provider_columns() -> None:
    specification = epcsaft_regression.FIGIEL_AQUEOUS_KIJ_V1
    models = _aqueous_kij_models(specification)
    capsules = tuple(native_sdk(model) for model in models)
    payload = _aqueous_kij_native_payload(specification)

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


def test_conditional_campaign_evidence_retains_separate_statuses() -> None:
    evidence_path = (
        Path(__file__).parents[1]
        / "evidence"
        / "figiel-aqueous-kij-conditional-recovery.json"
    )
    evidence_bytes = evidence_path.read_bytes()
    assert hashlib.sha256(evidence_bytes).hexdigest() == (
        "6796624b77152d04050293ae4ae6140bf54b0afa181999dc345d69dd7d8e5d51"
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
        "1a594a537e632d902878ab4d6c76ce04ac580ebb"
    )
    assert evidence["regression_artifact"]["wheel_sha256"] == (
        "a9283603156edb6f564af50b320ba0807ccee2e37e0feee42851d00ea278e221"
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
