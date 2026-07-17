from __future__ import annotations

import importlib.util
from pathlib import Path
from zipfile import ZipFile

import pytest

from epcsaft import EPCSAFT, ParameterBundle, native_sdk

import epcsaft_regression._native as native


def _provider_model() -> EPCSAFT:
    parameters = ParameterBundle.from_catalog(
        "gross-2001-methane-ethane", version=1
    ).select(("methane",))
    return EPCSAFT(parameters)


def test_parameterized_capsule_tail_is_validated_from_installed_provider() -> None:
    model = _provider_model()
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
    fake_wheel = tmp_path / "epcsaft_regression-0.1.0.dev0-py3-none-any.whl"
    with ZipFile(fake_wheel, "w") as wheel:
        wheel.writestr(
            "epcsaft_regression-0.1.0.dev0.dist-info/METADATA",
            "Metadata-Version: 2.1\nName: epcsaft-regression\nVersion: 0.1.0.dev0\n",
        )
        wheel.writestr("epcsaft_regression/__init__.py", "not the installed package\n")

    with pytest.raises(SystemExit, match="differs from wheel"):
        runner._require_installed_distribution_matches_wheel(
            fake_wheel,
            "epcsaft-regression",
        )
