from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
from zipfile import ZipFile

import pytest

from epcsaft import Mixture, Parameters, native_sdk

import epcsaft_regression._native as native


def _provider_model(component_id: str) -> Mixture:
    parameters = Parameters.from_catalog(
        "gross-2001-methane-ethane",
        components=(component_id,),
        version=1,
    )
    return Mixture(parameters)


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


def test_candidate_runners_bind_provider_frontend_0_2(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tools_path = Path(__file__).parents[1] / "tools"
    monkeypatch.syspath_prepend(str(tools_path))

    candidate_spec = importlib.util.spec_from_file_location(
        "run_candidate",
        tools_path / "run_candidate.py",
    )
    assert candidate_spec is not None and candidate_spec.loader is not None
    candidate = importlib.util.module_from_spec(candidate_spec)
    candidate_spec.loader.exec_module(candidate)

    born_spec = importlib.util.spec_from_file_location(
        "run_figiel_born_candidate",
        tools_path / "run_figiel_born_candidate.py",
    )
    assert born_spec is not None and born_spec.loader is not None
    born = importlib.util.module_from_spec(born_spec)
    born_spec.loader.exec_module(born)

    expected_commit = "06fb933e0b02ea87eb553a0a27d7a5ddb2077d72"
    expected_wheel = "1b8d69aba5f24936040de52eda6db1d7c8306b1cca38a48b105986fa6b657806"
    expected_header = "3c7bcce4748edb14b19c2ba486fd8f7ddddc0593ed333037fe06e09ae55237ba"

    assert candidate.PROVIDER_COMMIT == born.PROVIDER_COMMIT == expected_commit
    assert candidate.PROVIDER_WHEEL_SHA256 == born.PROVIDER_WHEEL_SHA256 == expected_wheel
    assert born.PROVIDER_HEADER_SHA256 == expected_header
    assert not hasattr(candidate, "PROVIDER_TEST_RECEIPT_SHA256")


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
