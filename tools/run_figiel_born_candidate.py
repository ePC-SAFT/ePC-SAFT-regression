from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path

import run_candidate as artifact_binding


PROVIDER_COMMIT = "907b077ec6f841a8a028fc759df14f899c79339c"
PROVIDER_TREE = "2b315113c9961a16f75c776783f704db54d75e44"
PROVIDER_WHEEL_SHA256 = "c327b9a176e54bfc79b625cca7f0c87f2a62fc7d87059826e40c9d70e214f0cd"
PROVIDER_HEADER_SHA256 = "610cc480f05c3e17e431d26fd1b2c8628eec3e2adb412102a284d4d5d6eb8171"
STAGED_PROVIDER_COMMIT = "238ff2f59b105126da059558958ca8b28749ad96"
STAGED_PROVIDER_TREE = "989d68822d22f52ab02ffedd7f26bcab6c78a1a3"
STAGED_PROVIDER_WHEEL_SHA256 = (
    "88b1b0ebda2212499cb0e270f8ee57e336c3c6ea833292cbcb57ecf2045cd029"
)
STAGED_PROVIDER_HEADER_SHA256 = (
    "5765af9fb4d90f70070eeeec12a8ccb63745f961d2b9abefbd8405220d291e67"
)


def _canonical_evidence_bytes(payload: dict[str, object]) -> bytes:
    record = dict(payload)
    record["evidence_payload_sha256"] = artifact_binding._canonical_json_sha256(payload)
    return (json.dumps(record, indent=2, sort_keys=True) + "\n").encode()


def _staged_conclusion(result: object, formula_checks: list[dict[str, object]]) -> str:
    formula_valid = bool(formula_checks) and all(
        check.get("passed") is True for check in formula_checks
    )
    if result.scientifically_valid and formula_valid:
        return "FIGIEL_STAGED_AQUEOUS_RECOVERY_CANDIDATE_PASSED"
    if (
        result.solver_converged
        and result.numerically_converged
        and result.physically_valid
        and result.workflow_valid
        and formula_valid
    ):
        return "SOURCE_DESCRIBED_STAGED_RECOVERY_DID_NOT_REPRODUCE_PRINTED_TUPLE"
    return "BLOCKED_STRICT_LOCAL_GATES"


def _run_staged_candidate(
    arguments: argparse.Namespace,
    provider_binding: dict[str, object],
    regression_binding: dict[str, object],
) -> int:
    import epcsaft
    import epcsaft_regression
    import epcsaft_regression._native as native
    from epcsaft import native_sdk
    from epcsaft_regression import (
        FIGIEL_STAGED_AQUEOUS_RECOVERY_V1,
        fit_figiel_staged_aqueous_parameters,
    )
    from epcsaft_regression.workflow import (
        AQUEOUS_KIJ_COLUMNS,
        _aqueous_native_payload,
        _figiel_models,
    )

    artifact_binding._require_import_origin(
        epcsaft.__file__, provider_binding, "epcsaft/__init__.py"
    )
    artifact_binding._require_import_origin(
        epcsaft_regression.__file__,
        regression_binding,
        "epcsaft_regression/__init__.py",
    )
    provider_header = (
        Path(epcsaft.__file__).parent / "include" / "epcsaft" / "native_sdk_v1.h"
    )
    artifact_binding._require_hash(
        provider_header, STAGED_PROVIDER_HEADER_SHA256, "provider header"
    )

    specification = FIGIEL_STAGED_AQUEOUS_RECOVERY_V1
    result = fit_figiel_staged_aqueous_parameters()
    _, salt_models = _figiel_models(
        result.born_diameters_angstrom,
        result.water_solvation_factor,
        result.aqueous_kij,
    )
    fingerprints = tuple(model.parameter_fingerprint for model in salt_models)
    capsules = tuple(native_sdk(model) for model in salt_models)
    formula_checks = []
    for stage, parameters in (
        ("solvation_factor", (result.water_solvation_factor,)),
        ("aqueous_kij", result.aqueous_kij),
    ):
        starts = (("evidence", tuple(parameters)),)
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
        matches = 0
        for row_index, row in enumerate(rows):
            ratio = float(row[5]) / float(row[3])
            columns = (
                (0,)
                if stage == "solvation_factor"
                else AQUEOUS_KIJ_COLUMNS[str(row[1])]
            )
            expected = [0.0] * parameter_count
            for column, derivative in zip(columns, row[7], strict=True):
                expected[column] = -ratio * float(derivative)
            observed = tuple(
                float(value)
                for value in jacobian[
                    row_index * parameter_count:(row_index + 1) * parameter_count
                ]
            )
            matches += (
                float(residuals[row_index]) == 1.0 - ratio
                and observed == tuple(expected)
            )
        formula_checks.append(
            {
                "stage": stage,
                "compiled_identity_round_trip": tuple(identity) == payload[0],
                "rows_checked": len(rows),
                "rows_with_exact_residual_and_jacobian_assembly": matches,
                "passed": matches == len(rows) and tuple(identity) == payload[0],
            }
        )

    conclusion = _staged_conclusion(result, formula_checks)
    subject = {
        "capability": specification.specification_id,
        "owner": "ePC-SAFT/ePC-SAFT-regression",
        "authority_status": "authority-neutral non-production package candidate",
        "conclusion": conclusion,
        "source": {
            "validation_commit": specification.source_validation_commit,
            "validation_tree": specification.source_validation_tree,
            "ledger_sha256": specification.source_ledger_sha256,
            "parameter_packet_sha256": specification.source_parameter_packet_sha256,
            "metadata_sha256": specification.source_metadata_sha256,
            "si_extraction_sha256": specification.source_si_extraction_sha256,
            "packaged_hamer_wu_sha256": specification.source_csv_sha256,
            "residual_row_count": 164,
            "stage_b_nabr_row_count": 21,
        },
        "provider": {
            "commit": STAGED_PROVIDER_COMMIT,
            "tree": STAGED_PROVIDER_TREE,
            "wheel": arguments.provider_wheel.name,
            "wheel_sha256": STAGED_PROVIDER_WHEEL_SHA256,
            "installed_header_sha256": STAGED_PROVIDER_HEADER_SHA256,
            "capsule": "epcsaft.native_sdk.v1",
            "derivative_order": 1,
        },
        "regression": {
            "commit": arguments.regression_commit,
            "tree": arguments.regression_tree,
            "wheel": arguments.regression_wheel.name,
            "wheel_sha256": artifact_binding._sha256(arguments.regression_wheel),
        },
        "installed_artifact_binding": {
            "method": "all non-RECORD wheel members matched installed files byte-for-byte",
            "provider": {
                key: value
                for key, value in provider_binding.items()
                if key != "verified_paths"
            },
            "regression": {
                key: value
                for key, value in regression_binding.items()
                if key != "verified_paths"
            },
        },
        "specification": asdict(specification),
        "exact_residual_jacobian_consumption": formula_checks,
        "result": asdict(result),
        "claim_limits": {
            "all_rows_are_training_data": True,
            "predictive_claim": False,
            "uncertainty_claim": False,
            "global_identifiability_claim": False,
            "provider_catalog_admission": False,
            "downstream_mea_readiness_claim": False,
        },
        "runner_checks": {
            "canonical_serialization_repeated_in_process": True,
        },
    }
    payload_record = {
        "subject": subject,
        "subject_sha256": artifact_binding._canonical_json_sha256(subject),
    }
    evidence = _canonical_evidence_bytes(payload_record)
    if evidence != _canonical_evidence_bytes(payload_record):
        raise RuntimeError("canonical staged evidence serialization was unstable")
    arguments.output.write_bytes(evidence)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate one Figiel candidate evidence record.")
    parser.add_argument("--workflow", choices=("born", "staged"), default="born")
    parser.add_argument("--provider-wheel", type=Path, required=True)
    parser.add_argument("--regression-wheel", type=Path, required=True)
    parser.add_argument("--regression-commit", required=True)
    parser.add_argument("--regression-tree", required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    expected_provider_wheel = (
        STAGED_PROVIDER_WHEEL_SHA256
        if arguments.workflow == "staged"
        else PROVIDER_WHEEL_SHA256
    )
    artifact_binding._require_hash(
        arguments.provider_wheel, expected_provider_wheel, "provider wheel"
    )
    provider_binding = artifact_binding._require_installed_distribution_matches_wheel(
        arguments.provider_wheel, "epcsaft"
    )
    regression_binding = artifact_binding._require_installed_distribution_matches_wheel(
        arguments.regression_wheel, "epcsaft-regression"
    )
    if arguments.workflow == "staged":
        return _run_staged_candidate(arguments, provider_binding, regression_binding)

    import epcsaft
    import epcsaft_regression
    import epcsaft_regression._native as native
    from epcsaft import EPCSAFT, ParameterBundle, native_sdk
    from epcsaft_regression import (
        FIGIEL_BORN_DIAMETER_TRACER_V1,
        fit_figiel_born_diameters,
    )
    from epcsaft_regression.workflow import _born_native_payload

    artifact_binding._require_import_origin(
        epcsaft.__file__, provider_binding, "epcsaft/__init__.py"
    )
    artifact_binding._require_import_origin(
        epcsaft_regression.__file__,
        regression_binding,
        "epcsaft_regression/__init__.py",
    )
    provider_header = (
        Path(epcsaft.__file__).parent / "include" / "epcsaft" / "native_sdk_v1.h"
    )
    artifact_binding._require_hash(provider_header, PROVIDER_HEADER_SHA256, "provider header")

    specification = FIGIEL_BORN_DIAMETER_TRACER_V1
    catalog = ParameterBundle.from_catalog("figiel-2025-reference-electrolytes", version=1)
    models = tuple(
        EPCSAFT(catalog.select(target.component_order)) for target in specification.targets
    )
    capsules = tuple(native_sdk(model) for model in models)
    payload = _born_native_payload(specification)
    trial = (2.9, 3.2, 4.3, 3.9, 4.7)
    trial_rows = native.evaluate_born(capsules, payload, trial)[2]
    derivative_checks = []
    for index, target in enumerate(specification.targets):
        differences = []
        for step in (1.0e-4, 5.0e-5):
            plus = list(trial)
            minus = list(trial)
            plus[index] += step
            minus[index] -= step
            value_plus = native.evaluate_born(capsules, payload, tuple(plus))[2][index][0]
            value_minus = native.evaluate_born(capsules, payload, tuple(minus))[2][index][0]
            differences.append((value_plus - value_minus) / (2.0 * step))
        exact = float(trial_rows[index][1])
        tolerance = max(
            1.0e-3,
            20.0 * abs(differences[0] - differences[1]),
            2.0e-8 * abs(exact),
        )
        derivative_checks.append(
            {
                "target_id": target.target_id,
                "trial_diameter_angstrom": trial[index],
                "exact_j_per_mol_per_angstrom": exact,
                "centered_h_j_per_mol_per_angstrom": differences[0],
                "centered_h_over_2_j_per_mol_per_angstrom": differences[1],
                "absolute_difference": abs(exact - differences[1]),
                "acceptance_tolerance": tolerance,
                "passed": abs(exact - differences[1]) <= tolerance,
            }
        )

    result = fit_figiel_born_diameters(models=models)
    conclusion = (
        "FIGIEL_BORN_TRACER_CANDIDATE_PASSED"
        if result.scientifically_valid
        else "BLOCKED_STRICT_LOCAL_GATES"
    )
    subject = {
        "capability": "figiel-2025-five-ion-born-diameter-tracer-v1",
        "owner": "ePC-SAFT/ePC-SAFT-regression",
        "authority_status": "authority-neutral non-production package candidate",
        "conclusion": conclusion,
        "strict_local_gates_passed": result.scientifically_valid,
        "source": {
            "validation_commit": specification.source_validation_commit,
            "validation_tree": specification.source_validation_tree,
            "ledger_sha256": specification.source_ledger_sha256,
            "parameter_packet_sha256": specification.source_parameter_packet_sha256,
            "metadata_sha256": specification.source_metadata_sha256,
            "packaged_five_target_sha256": specification.packaged_targets_sha256,
            "residual_target_count": 5,
            "underlying_support_rows_copied_or_fitted": 0,
        },
        "provider": {
            "commit": PROVIDER_COMMIT,
            "tree": PROVIDER_TREE,
            "wheel": arguments.provider_wheel.name,
            "wheel_sha256": PROVIDER_WHEEL_SHA256,
            "installed_header_sha256": PROVIDER_HEADER_SHA256,
            "capsule": "epcsaft.native_sdk.v1",
            "entry": "evaluate_ion_solvation_born",
            "derivative_order": 1,
        },
        "regression": {
            "commit": arguments.regression_commit,
            "tree": arguments.regression_tree,
            "wheel": arguments.regression_wheel.name,
            "wheel_sha256": artifact_binding._sha256(arguments.regression_wheel),
        },
        "installed_artifact_binding": {
            "method": "all non-RECORD wheel members matched installed files byte-for-byte",
            "provider": {
                key: value for key, value in provider_binding.items() if key != "verified_paths"
            },
            "regression": {
                key: value for key, value in regression_binding.items() if key != "verified_paths"
            },
        },
        "specification": asdict(specification),
        "derivative_checks": derivative_checks,
        "result": asdict(result),
        "claim_limits": {
            "all_targets_are_training_equations": True,
            "predictive_claim": False,
            "uncertainty_claim": False,
            "global_identifiability_claim": False,
            "provider_catalog_admission": False,
        },
    }
    payload_record = {
        "subject": subject,
        "subject_sha256": artifact_binding._canonical_json_sha256(subject),
    }
    arguments.output.write_bytes(_canonical_evidence_bytes(payload_record))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
