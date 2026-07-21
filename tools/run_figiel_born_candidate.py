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


def _canonical_evidence_bytes(payload: dict[str, object]) -> bytes:
    record = dict(payload)
    record["evidence_payload_sha256"] = artifact_binding._canonical_json_sha256(payload)
    return (json.dumps(record, indent=2, sort_keys=True) + "\n").encode()


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate the Figiel Born-tracer candidate evidence.")
    parser.add_argument("--provider-wheel", type=Path, required=True)
    parser.add_argument("--regression-wheel", type=Path, required=True)
    parser.add_argument("--regression-commit", required=True)
    parser.add_argument("--regression-tree", required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    artifact_binding._require_hash(
        arguments.provider_wheel, PROVIDER_WHEEL_SHA256, "provider wheel"
    )
    provider_binding = artifact_binding._require_installed_distribution_matches_wheel(
        arguments.provider_wheel, "epcsaft"
    )
    regression_binding = artifact_binding._require_installed_distribution_matches_wheel(
        arguments.regression_wheel, "epcsaft-regression"
    )

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
        else "BLOCKED_PUBLISHED_DIAMETER_RECOVERY"
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
