from __future__ import annotations

import argparse
from dataclasses import asdict
from email.parser import BytesParser
import hashlib
from importlib import metadata
import json
from pathlib import Path
import platform
from zipfile import ZipFile


PROVIDER_COMMIT = "4b10cb899c94687cae734980285badb224dc95e6"
PROVIDER_WHEEL_SHA256 = "f92f79c8d6f614660e5c201b7061c9b02b5cd1a25a4ed8c8fee0b59adaabf2bf"
PROVIDER_TEST_RECEIPT_SHA256 = "07447721abaca946c6e9221e7d49e431e13fcb8e6867944f67b6ba8337901480"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _require_hash(path: Path, expected: str, label: str) -> None:
    actual = _sha256(path)
    if actual != expected:
        raise SystemExit(f"{label} SHA-256 mismatch: expected {expected}, got {actual}")


def _normalized_distribution_name(name: str) -> str:
    return name.lower().replace("-", "_").replace(".", "_")


def _require_installed_distribution_matches_wheel(
    wheel_path: Path,
    distribution_name: str,
) -> dict[str, object]:
    """Prove every wheel member used at runtime is installed byte-for-byte."""
    distribution = metadata.distribution(distribution_name)
    normalized_expected_name = _normalized_distribution_name(distribution_name)
    verified_paths: dict[str, str] = {}
    with ZipFile(wheel_path) as wheel:
        members = [item for item in wheel.infolist() if not item.is_dir()]
        if not members:
            raise SystemExit(f"{wheel_path} contains no files")
        for member in members:
            member_path = Path(member.filename)
            if member_path.is_absolute() or ".." in member_path.parts:
                raise SystemExit(f"unsafe wheel member path: {member.filename}")
        metadata_members = [
            item
            for item in members
            if item.filename.endswith(".dist-info/METADATA")
        ]
        if len(metadata_members) != 1:
            raise SystemExit(f"{wheel_path} must contain exactly one METADATA file")
        wheel_metadata = BytesParser().parsebytes(wheel.read(metadata_members[0]))
        wheel_name = wheel_metadata["Name"]
        wheel_version = wheel_metadata["Version"]
        if _normalized_distribution_name(wheel_name) != normalized_expected_name:
            raise SystemExit(
                f"wheel distribution mismatch: expected {distribution_name}, got {wheel_name}"
            )
        if distribution.version != wheel_version:
            raise SystemExit(
                f"installed {distribution_name} version {distribution.version} "
                f"does not match wheel version {wheel_version}"
            )
        for member in members:
            if member.filename.endswith(".dist-info/RECORD"):
                continue
            installed_path = Path(distribution.locate_file(member.filename)).resolve()
            if not installed_path.is_file():
                raise SystemExit(
                    f"installed {distribution_name} is missing wheel member {member.filename}"
                )
            wheel_sha256 = hashlib.sha256(wheel.read(member)).hexdigest()
            installed_sha256 = _sha256(installed_path)
            if installed_sha256 != wheel_sha256:
                raise SystemExit(
                    f"installed {distribution_name} member differs from wheel: "
                    f"{member.filename}"
                )
            verified_paths[member.filename] = str(installed_path)
    return {
        "distribution": distribution.metadata["Name"],
        "version": distribution.version,
        "verified_member_count": len(verified_paths),
        "verified_paths": verified_paths,
    }


def _require_import_origin(
    module_file: str,
    binding: dict[str, object],
    expected_member: str,
) -> None:
    verified_paths = binding["verified_paths"]
    if not isinstance(verified_paths, dict):
        raise SystemExit("invalid installed-distribution binding evidence")
    expected_path = verified_paths.get(expected_member)
    if expected_path is None or Path(module_file).resolve() != Path(expected_path):
        raise SystemExit(
            f"import origin {module_file} is not the verified wheel member {expected_member}"
        )


def _rms(values: list[float]) -> float:
    return (sum(value * value for value in values) / len(values)) ** 0.5


def _error_metrics(rows: list[object], interpretation: str) -> dict[str, object]:
    valid_rows = [row for row in rows if row.physically_valid]
    pressure_errors = [row.pressure_relative_error for row in valid_rows]
    density_errors = [row.liquid_density_relative_error for row in valid_rows]
    return {
        "requested_temperatures_k": [row.temperature_k for row in rows],
        "valid_temperatures_k": [row.temperature_k for row in valid_rows],
        "failed_rows": [
            {
                "temperature_k": row.temperature_k,
                "failure_reasons": list(row.failure_reasons),
            }
            for row in rows
            if not row.physically_valid
        ],
        "pressure_relative_error_rms": _rms(pressure_errors) if valid_rows else None,
        "pressure_relative_error_max_abs": max(map(abs, pressure_errors)) if valid_rows else None,
        "liquid_density_relative_error_rms": _rms(density_errors) if valid_rows else None,
        "liquid_density_relative_error_max_abs": max(map(abs, density_errors))
        if valid_rows
        else None,
        "acceptance_threshold": None,
        "interpretation": interpretation,
    }


def _canonical_json_sha256(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _canonical_receipt_bytes(payload: dict[str, object]) -> bytes:
    receipt = dict(payload)
    receipt["receipt_payload_sha256"] = _canonical_json_sha256(payload)
    return (json.dumps(receipt, indent=2, sort_keys=True) + "\n").encode("utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate one pure-saturation regression receipt.")
    parser.add_argument("--component", choices=("methane", "ethane"), required=True)
    parser.add_argument("--provider-wheel", type=Path, required=True)
    parser.add_argument("--provider-test-receipt", type=Path, required=True)
    parser.add_argument("--regression-wheel", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--reviewer")
    parser.add_argument("--review-path")
    parser.add_argument("--review-sha256")
    arguments = parser.parse_args()
    _require_hash(arguments.provider_wheel, PROVIDER_WHEEL_SHA256, "provider wheel")
    _require_hash(
        arguments.provider_test_receipt,
        PROVIDER_TEST_RECEIPT_SHA256,
        "provider test receipt",
    )
    regression_wheel_sha256 = _sha256(arguments.regression_wheel)
    repository_root = Path(__file__).resolve().parents[1]
    review_arguments = (arguments.reviewer, arguments.review_path, arguments.review_sha256)
    if any(review_arguments) and not all(review_arguments):
        raise SystemExit("reviewer, review-path, and review-sha256 must be supplied together")
    review_record: dict[str, object]
    if all(review_arguments):
        review_path = Path(arguments.review_path).resolve()
        try:
            portable_review_path = review_path.relative_to(repository_root).as_posix()
        except ValueError as error:
            raise SystemExit("independent review must be inside the regression repository") from error
        _require_hash(review_path, arguments.review_sha256, "independent review")
        review_record = {
            "status": "provided",
            "identity": arguments.reviewer,
            "path": portable_review_path,
            "sha256": arguments.review_sha256,
        }
    else:
        review_record = {"status": "pending"}
    provider_binding = _require_installed_distribution_matches_wheel(
        arguments.provider_wheel,
        "epcsaft",
    )
    regression_binding = _require_installed_distribution_matches_wheel(
        arguments.regression_wheel,
        "epcsaft-regression",
    )

    import epcsaft
    import epcsaft_regression
    from epcsaft import EPCSAFT, ParameterBundle
    from epcsaft_regression import (
        ETHANE_SATURATION_FIT_V1,
        METHANE_SATURATION_FIT_V1,
        fit_pure_saturation,
        load_pure_saturation_dataset,
    )

    _require_import_origin(
        epcsaft.__file__,
        provider_binding,
        "epcsaft/__init__.py",
    )
    _require_import_origin(
        epcsaft_regression.__file__,
        regression_binding,
        "epcsaft_regression/__init__.py",
    )

    component_id = arguments.component
    dataset = load_pure_saturation_dataset(component_id)
    specification = (
        METHANE_SATURATION_FIT_V1
        if component_id == "methane"
        else ETHANE_SATURATION_FIT_V1
    )
    parameters = ParameterBundle.from_catalog(
        "gross-2001-methane-ethane", version=1
    ).select((component_id,))
    model = EPCSAFT(parameters)
    result = fit_pure_saturation(
        model=model,
        dataset=dataset,
        specification=specification,
    )
    if not (result.solver_converged and result.numerically_converged and result.physically_valid):
        raise SystemExit(f"candidate failed strict gates: {result.failure_reasons}")

    held_out_temperatures = frozenset(dataset.held_out_temperatures_k)
    stress_temperatures = frozenset(dataset.stress_temperatures_k)
    training = [row for row in result.reporting_rows if row.training]
    held_out = [row for row in result.reporting_rows if row.temperature_k in held_out_temperatures]
    stress = [row for row in result.reporting_rows if row.temperature_k in stress_temperatures]
    problem = asdict(specification)
    source = asdict(dataset.source)
    rows = [asdict(row) for row in dataset.rows]
    result_record = asdict(result)
    subject = {
        "capability": f"pure-{component_id}-saturation-parameter-candidate-v1",
        "component_id": component_id,
        "owner": "ePC-SAFT/ePC-SAFT-regression",
        "authority_status": (
            "promotion-0020 accepted for the exact reproducible workflow; replay changes no authority"
            if component_id == "methane"
            else "local candidate; independent review and validation admission pending"
        ),
        "source": source,
        "rows": rows,
        "training_row_ids": list(dataset.training_row_ids),
        "held_out_row_ids": [row.row_id for row in dataset.held_out_rows],
        "stress_row_ids": [row.row_id for row in dataset.stress_rows],
        "provider": {
            "commit": PROVIDER_COMMIT,
            "wheel": arguments.provider_wheel.name,
            "wheel_sha256": PROVIDER_WHEEL_SHA256,
            "test_receipt": arguments.provider_test_receipt.name,
            "test_receipt_sha256": PROVIDER_TEST_RECEIPT_SHA256,
            "capsule": "epcsaft.native_sdk.v1",
            "entry": "evaluate_pure_phase_parameters",
            "coordinate_order": [
                "amount_mol",
                "volume_m3",
                "segment_count",
                "segment_diameter_angstrom",
                "dispersion_energy_over_k_kelvin",
            ],
            "source_parameter_fingerprint": model.parameter_fingerprint,
        },
        "regression_artifact": {
            "wheel": arguments.regression_wheel.name,
            "wheel_sha256": regression_wheel_sha256,
        },
        "installed_artifact_binding": {
            "method": "all non-RECORD wheel members matched installed files byte-for-byte before package import",
            "provider": {
                key: value for key, value in provider_binding.items() if key != "verified_paths"
            },
            "regression": {
                key: value for key, value in regression_binding.items() if key != "verified_paths"
            },
            "import_origins": {
                "provider": "epcsaft/__init__.py",
                "regression": "epcsaft_regression/__init__.py",
            },
        },
        "problem": problem,
        "transforms": {
            "parameters": "p_j = start_j + parameter_scale_j * z_j",
            "row_volumes": "V_phase = V_phase_start * exp(u_phase)",
            "reporting_pressure": "P_report = P_observed * exp(u_pressure)",
            "volume_starts": {
                "liquid": "molar_mass / observed_liquid_mass_density",
                "vapor": "R*T / observed_pressure",
            },
            "gas_constant_j_per_mol_k": 8.31446261815324,
            "start_provenance": (
                "retained M5 displaced methane start: (1.0*1.08, 3.7039*0.96, 150.03*1.05)"
                if component_id == "methane"
                else "Gross 2001 ethane source parameters: (1.6069, 3.5206, 191.42)"
            ),
            "molar_mass_provenance": (
                "fixed formulation constant 0.016043 kg/mol retained by the admitted M5 plan"
                if component_id == "methane"
                else "source-backed ethane molar mass 0.030070 kg/mol"
            ),
        },
        "residuals": {
            "order": list(specification.residual_names),
            "raw": [
                "P_liquid - P_observed [Pa]",
                "P_vapor - P_observed [Pa]",
                "mu_liquid/(RT) - mu_vapor/(RT) [1]",
                "molar_mass/V_liquid - observed_liquid_mass_density [kg/m3]",
            ],
            "scales": ["P_observed", "P_observed", 1.0, "rho_liquid_observed"],
            "weights": list(specification.residual_weights),
            "weight_meaning": "equal per-row normalization; not an uncertainty claim",
            "jacobian": "provider Hessian plus exact affine/log/linear chain rules",
        },
        "result": result_record,
        "training_reporting": _error_metrics(
            training,
            "fit-row descriptive errors; objective residual structure remains the solver evidence",
        ),
        "held_out_reporting": _error_metrics(
            held_out,
            "descriptive candidate evidence; validation owns admission",
        ),
        "stress_reporting": _error_metrics(
            stress,
            "reported domain stress only; cannot determine solver or scientific acceptance",
        ),
        "execution": {
            "artifacts": {
                "provider_wheel": arguments.provider_wheel.name,
                "provider_test_receipt": arguments.provider_test_receipt.name,
                "regression_wheel": arguments.regression_wheel.name,
            },
            "python": platform.python_version(),
            "platform": platform.platform(),
            "ceres": "2.2.0 system package used by the linked native target",
        },
        "exclusions": [
            "binary k_ij",
            "association parameters",
            "electrolyte targets",
            "reactive targets",
            "generic regression families",
            "provider catalog persistence",
            "global identifiability claim",
            "predictive validation admission",
            "runtime authority promotion",
        ],
        "evidence_boundary": {
            "reporting_block_directional_jacobian": (
                "promotion-receipt evidence item; the existing private native test surface "
                "does not expose reporting Jacobians and no production seam was added"
            ),
        },
    }
    subject_sha256 = _canonical_json_sha256(subject)
    receipt = {
        "schema_version": "epcsaft-regression-candidate-fit-receipt-v1",
        "status": "accepted-workflow-replay" if component_id == "methane" else "candidate",
        "subject_sha256": subject_sha256,
        "subject": subject,
        "producer": "codex:/root",
        "independent_reviewer": review_record,
        "approval": {
            "implementation_scope": "approved in delegated user prompt",
            "runtime_authority_change": "not requested; not granted",
            "validation_admission": "pending",
        },
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_bytes(_canonical_receipt_bytes(receipt))
    print(json.dumps({"output": str(arguments.output), "subject_sha256": subject_sha256}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
