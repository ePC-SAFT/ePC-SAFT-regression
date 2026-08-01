from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
from zipfile import ZipFile

import epcsaft_regression._native as native
import pytest
from epcsaft import Mixture, Parameters, native_sdk
from parameter_cases import neutral_parameters


def _provider_model(component_id: str) -> Mixture:
    return Mixture(Parameters.from_dictionary(neutral_parameters((component_id,))))


@pytest.mark.parametrize("component_id", ("methane", "ethane"))
def test_parameterized_capsule_tail_is_validated_from_installed_provider(
    component_id: str,
) -> None:
    model = _provider_model(component_id)
    capsule = native_sdk(model)

    abi_version, table_size, parameterized_result_size, fingerprint = native.transport_info(
        capsule
    )

    assert abi_version == 1
    assert table_size >= native.minimum_parameterized_table_size()
    assert parameterized_result_size == native.parameterized_result_size()
    assert fingerprint == model.parameter_fingerprint


def test_provider_header_is_from_installed_wheel_not_a_sibling_source_tree() -> None:
    import epcsaft

    header = Path(epcsaft.__file__).parent / "include" / "epcsaft" / "native_sdk_v1.h"

    assert header.is_file()
    assert "site-packages" in header.as_posix()
    assert "/ePC-SAFT-project/ePC-SAFT/src/" not in header.as_posix()


def test_receipt_runner_rejects_wheel_that_differs_from_installed_runtime(
    tmp_path: Path,
) -> None:
    runner_path = Path(__file__).parents[1] / "tools" / "run_candidate.py"
    module_spec = importlib.util.spec_from_file_location("candidate_runner_test", runner_path)
    assert module_spec is not None and module_spec.loader is not None
    runner = importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(runner)
    fake_wheel = tmp_path / "epcsaft_regression-0.2.0.dev0-py3-none-any.whl"
    with ZipFile(fake_wheel, "w") as wheel:
        wheel.writestr(
            "epcsaft_regression-0.2.0.dev0.dist-info/METADATA",
            "Metadata-Version: 2.1\nName: epcsaft-regression\nVersion: 0.2.0.dev0\n",
        )
        wheel.writestr("epcsaft_regression/__init__.py", "not the installed package\n")

    with pytest.raises(SystemExit, match="differs from wheel"):
        runner._require_installed_distribution_matches_wheel(
            fake_wheel,
            "epcsaft-regression",
        )


def test_candidate_runners_share_final_provider_identity(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tools_path = Path(__file__).parents[1] / "tools"
    monkeypatch.syspath_prepend(str(tools_path))

    candidate_spec = importlib.util.spec_from_file_location(
        "run_candidate_identity", tools_path / "run_candidate.py"
    )
    assert candidate_spec is not None and candidate_spec.loader is not None
    candidate = importlib.util.module_from_spec(candidate_spec)
    candidate_spec.loader.exec_module(candidate)

    born_spec = importlib.util.spec_from_file_location(
        "run_figiel_born_identity", tools_path / "run_figiel_born_candidate.py"
    )
    assert born_spec is not None and born_spec.loader is not None
    born = importlib.util.module_from_spec(born_spec)
    born_spec.loader.exec_module(born)

    expected_commit = "7b97bab039e1c50a6f89522698af80493bea5f9e"
    expected_tree = "d082a8f102b32705b6cd6669a3e31a8d4ea8acd0"
    expected_wheel = "1567cda72e1b525526dc0e647af0c6fe711edcb70bc4cee08f06284e847956d9"
    expected_header = "881f5ec87293de8b1f3c25c16018aa94be69775fede2ec5426fcbb08e257fecd"
    expected_library = "fd624add206b8d783cd079db320b6dba64083063af2f29faf5ce82d1cf4743eb"

    assert candidate.PROVIDER_COMMIT == born.PROVIDER_COMMIT == expected_commit
    assert candidate.PROVIDER_TREE == born.PROVIDER_TREE == expected_tree
    assert candidate.PROVIDER_WHEEL_SHA256 == born.PROVIDER_WHEEL_SHA256 == expected_wheel
    assert candidate.PROVIDER_HEADER_SHA256 == born.PROVIDER_HEADER_SHA256 == expected_header
    assert candidate.PROVIDER_LIBRARY_SHA256 == born.PROVIDER_LIBRARY_SHA256 == expected_library


def test_parameter_bundle_tree_hash_binds_paths_and_bytes(tmp_path: Path) -> None:
    runner_path = Path(__file__).parents[1] / "tools" / "run_candidate.py"
    module_spec = importlib.util.spec_from_file_location(
        "candidate_runner_bundle_hash", runner_path
    )
    assert module_spec is not None and module_spec.loader is not None
    runner = importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(runner)
    bundle = tmp_path / "parameters"
    nested = bundle / "nested"
    nested.mkdir(parents=True)
    (bundle / "bundle.toml").write_text('schema = "epcsaft.parameters"\n')
    (nested / "single.csv").write_text("record_id,value\nmethane,1\n")

    original = runner._directory_sha256(bundle)
    assert original == runner._directory_sha256(bundle)

    (nested / "single.csv").write_text("record_id,value\nmethane,2\n")
    assert runner._directory_sha256(bundle) != original


def test_candidate_receipt_has_one_canonical_reproducible_subject() -> None:
    runner_path = Path(__file__).parents[1] / "tools" / "run_candidate.py"
    module_spec = importlib.util.spec_from_file_location("candidate_runner_canonical", runner_path)
    assert module_spec is not None and module_spec.loader is not None
    runner = importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(runner)
    receipt_path = Path(__file__).parents[1] / "evidence" / "candidate-fit-receipt.json"
    receipt_bytes = receipt_path.read_bytes()
    receipt = json.loads(receipt_bytes)
    payload = dict(receipt)
    payload.pop("receipt_payload_sha256")

    assert runner._canonical_receipt_bytes(payload) == receipt_bytes
    assert runner._canonical_json_sha256(receipt["subject"]) == receipt["subject_sha256"]
    assert receipt["independent_reviewer"]["path"] == "docs/reviews/independent-review.md"
    assert not {"source", "rows", "training_row_ids", "problem"}.intersection(
        receipt.keys()
    )
    assert {"source", "rows", "training_row_ids", "problem"} <= receipt["subject"].keys()


def test_figiel_candidate_evidence_has_one_canonical_passed_subject() -> None:
    evidence_path = (
        Path(__file__).parents[1] / "evidence" / "figiel-born-diameter-candidate.json"
    )
    evidence_bytes = evidence_path.read_bytes()
    evidence = json.loads(evidence_bytes)
    payload = dict(evidence)
    payload.pop("evidence_payload_sha256")
    canonical_payload = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    expected_bytes = dict(payload)
    expected_bytes["evidence_payload_sha256"] = hashlib.sha256(canonical_payload).hexdigest()

    assert (json.dumps(expected_bytes, indent=2, sort_keys=True) + "\n").encode() == evidence_bytes
    subject = evidence["subject"]
    canonical_subject = json.dumps(
        subject, sort_keys=True, separators=(",", ":")
    ).encode()
    assert evidence["subject_sha256"] == hashlib.sha256(canonical_subject).hexdigest()
    assert subject["conclusion"] == "FIGIEL_BORN_TRACER_CANDIDATE_PASSED"
    assert subject["source"]["residual_target_count"] == 5
    assert subject["source"]["underlying_support_rows_copied_or_fitted"] == 0
    assert all(check["passed"] for check in subject["derivative_checks"])
    assert subject["result"]["solver_converged"]
    assert subject["result"]["numerically_converged"]
    assert subject["result"]["workflow_valid"]
    assert subject["result"]["scientifically_valid"]
    assert max(
        abs(parameter["published_delta_angstrom"])
        for parameter in subject["result"]["parameters"]
    ) > 0.0005
