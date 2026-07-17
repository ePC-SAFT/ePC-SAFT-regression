from __future__ import annotations

import ctypes
import math
from types import SimpleNamespace

import pytest
from epcsaft import EPCSAFT, ParameterBundle, native_sdk

import epcsaft_regression._native as native
from epcsaft_regression.records import (
    ETHANE_SATURATION_FIT_V1,
    METHANE_SATURATION_FIT_V1,
    PureSaturationFitSpecification,
    load_pure_saturation_dataset,
)
import epcsaft_regression.workflow as workflow
from epcsaft_regression.workflow import _native_payload, fit_pure_saturation


SPECIFICATIONS = {
    "methane": METHANE_SATURATION_FIT_V1,
    "ethane": ETHANE_SATURATION_FIT_V1,
}


class _ParameterizedResult(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("status", ctypes.c_int32),
        ("helmholtz", ctypes.c_double),
        ("gradient", ctypes.c_double * 5),
        ("hessian", ctypes.c_double * 25),
        ("pressure", ctypes.c_double),
        ("chemical_potential", ctypes.c_double),
        ("fingerprint", ctypes.c_char * 72),
        ("error", ctypes.c_char * 160),
    )


_ParameterizedCallback = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(_ParameterizedResult),
)


class _NativeSdkTable(ctypes.Structure):
    _fields_ = (
        ("abi_version", ctypes.c_uint32),
        ("table_size", ctypes.c_size_t),
        ("result_size", ctypes.c_size_t),
        ("model_context", ctypes.c_void_p),
        ("evaluate_pure_phase", ctypes.c_void_p),
        ("parameterized_result_size", ctypes.c_size_t),
        ("evaluate_pure_phase_parameters", ctypes.c_void_p),
    )


def _model(component_id: str) -> EPCSAFT:
    parameters = ParameterBundle.from_catalog(
        "gross-2001-methane-ethane", version=1
    ).select((component_id,))
    return EPCSAFT(parameters)


def _capsule(component_id: str) -> object:
    return native_sdk(_model(component_id))


def _payload(component_id: str) -> tuple[object, ...]:
    dataset = load_pure_saturation_dataset(component_id)
    specification = SPECIFICATIONS[component_id]
    model = _model(component_id)
    return _native_payload(dataset, specification, model.parameter_fingerprint)


def _failing_provider_capsule() -> tuple[object, _NativeSdkTable, object]:
    @_ParameterizedCallback
    def fail_evaluation(
        _context: object,
        _temperature: float,
        _amount: float,
        _volume: float,
        _segment_count: float,
        _segment_diameter: float,
        _dispersion_energy: float,
        result: ctypes.POINTER(_ParameterizedResult),
    ) -> int:
        result.contents.status = 3
        result.contents.error = b"synthetic provider domain failure"
        return 3

    table = _NativeSdkTable(
        abi_version=1,
        table_size=ctypes.sizeof(_NativeSdkTable),
        result_size=0,
        model_context=1,
        evaluate_pure_phase=None,
        parameterized_result_size=ctypes.sizeof(_ParameterizedResult),
        evaluate_pure_phase_parameters=ctypes.cast(fail_evaluation, ctypes.c_void_p).value,
    )
    capsule_new = ctypes.pythonapi.PyCapsule_New
    capsule_new.argtypes = (ctypes.c_void_p, ctypes.c_char_p, ctypes.c_void_p)
    capsule_new.restype = ctypes.py_object
    capsule = capsule_new(ctypes.addressof(table), b"epcsaft.native_sdk.v1", None)
    return capsule, table, fail_evaluation


def test_methane_start_residuals_match_accepted_provider_anchor() -> None:
    residuals, jacobian, diagnostics, fingerprint = native.evaluate(
        _capsule("methane"), _payload("methane"), (0.0,) * 11
    )

    assert residuals[4:8] == pytest.approx(
        (-30.793189605316272, -0.043171211223460487, -0.64887647749645083, 0.0),
        rel=2.0e-13,
        abs=2.0e-13,
    )
    assert len(residuals) == 16
    assert len(jacobian) == 16 * 11
    assert len(diagnostics) == 4
    assert fingerprint.startswith("sha256:")


@pytest.mark.parametrize("component_id", ("methane", "ethane"))
def test_exact_jacobian_matches_independent_directional_residual_difference(
    component_id: str,
) -> None:
    capsule = _capsule(component_id)
    payload = _payload(component_id)
    variables = (0.0,) * 11
    direction = (0.2, -0.1, 0.05, 0.01, -0.02, -0.015, 0.012, 0.008, -0.01, -0.006, 0.014)
    residuals, jacobian, _, _ = native.evaluate(capsule, payload, variables)
    step = 1.0e-6
    plus = tuple(value + step * delta for value, delta in zip(variables, direction, strict=True))
    minus = tuple(value - step * delta for value, delta in zip(variables, direction, strict=True))
    residuals_plus = native.evaluate(capsule, payload, plus)[0]
    residuals_minus = native.evaluate(capsule, payload, minus)[0]
    finite_difference = tuple(
        (right - left) / (2.0 * step)
        for right, left in zip(residuals_plus, residuals_minus, strict=True)
    )
    product = tuple(
        math.fsum(jacobian[row * 11 + column] * direction[column] for column in range(11))
        for row in range(16)
    )

    assert residuals
    assert product == pytest.approx(finite_difference, rel=2.0e-6, abs=2.0e-7)


@pytest.mark.parametrize("component_id", ("methane", "ethane"))
def test_public_workflow_returns_strict_component_diagnostics(
    component_id: str,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    model = _model(component_id)
    dataset = load_pure_saturation_dataset(component_id)
    specification = SPECIFICATIONS[component_id]
    native_solve = native.solve
    captured_transport: list[tuple[object, ...]] = []

    def record_transport(*args: object) -> tuple[object, ...]:
        transport = native_solve(*args)
        captured_transport.append(transport)
        return transport

    monkeypatch.setattr(native, "solve", record_transport)
    result = fit_pure_saturation(
        model=model,
        dataset=dataset,
        specification=specification,
    )

    assert result.component_id == component_id
    assert result.solver_converged, result.failure_reasons
    assert result.numerically_converged, result.failure_reasons
    assert result.physically_valid, result.failure_reasons
    assert result.termination == "CONVERGENCE"
    assert result.solution_usable
    assert result.final_cost < result.initial_cost
    assert result.provider_fingerprint == model.parameter_fingerprint
    assert result.provider_fingerprint in result.compiled_problem_identity
    assert len(result.parameters) == 3
    assert any(abs(item.movement) > 1.0e-8 for item in result.parameters)
    assert not any(item.active_bound for item in result.parameters)
    assert result.jacobian.complete_columns
    assert result.jacobian.full_rank == 11
    assert result.jacobian.parameter_rank == 3
    assert len(result.training_rows) == 4
    assert len(result.reporting_rows) == len(dataset.rows)
    assert tuple(row.temperature_k for row in result.reporting_rows) == tuple(
        row.temperature_k for row in dataset.rows
    )
    assert sum(row.training for row in result.reporting_rows) == 4
    assert all(
        row.physically_valid for row in result.reporting_rows if row.partition != "stress"
    )
    assert tuple(row.partition for row in result.reporting_rows) == tuple(
        "training"
        if row.temperature_k in dataset.training_temperatures_k
        else "held_out"
        if row.temperature_k in dataset.held_out_temperatures_k
        else "stress"
        for row in dataset.rows
    )
    assert result.confirmation_parameter_scaled_max_delta <= 1.0e-5
    assert result.confirmation_cost_relative_delta <= 1.0e-8
    assert result.failure_reasons == ()
    assert len(captured_transport) == 1
    assert len(captured_transport[0]) == 24
    assert tuple(captured_transport[0][22]) == result.compiled_problem_identity
    assert captured_transport[0][23] == ""


def test_generalized_workflow_preserves_accepted_methane_numerical_result() -> None:
    result = fit_pure_saturation(
        model=_model("methane"),
        dataset=load_pure_saturation_dataset("methane"),
        specification=METHANE_SATURATION_FIT_V1,
    )
    expected_reporting = (
        (100.0, 34626.07915160773, 436.84483289550474),
        (110.0, 88224.60866583801, 423.3969791365449),
        (120.0, 191083.41254773695, 409.3407614313209),
        (130.0, 366384.5067925305, 394.34823157350235),
        (140.0, 639981.9267634666, 377.95195367291694),
        (150.0, 1039603.4624909018, 359.43578939218077),
        (160.0, 1594405.452356535, 337.5806133970321),
        (170.0, 2334648.4404434026, 309.96652881184707),
        (180.0, 3290375.174877589, 270.4239126820564),
    )

    assert tuple(item.final for item in result.parameters) == pytest.approx(
        (0.9932081279826167, 3.717121437945618, 150.4888402511307),
        rel=2.0e-11,
        abs=2.0e-11,
    )
    assert result.initial_cost == pytest.approx(14340.021563034428, rel=2.0e-12)
    assert result.final_cost == pytest.approx(4.798586497669576e-6, rel=2.0e-9)
    assert result.jacobian.full_rank == 11
    assert result.jacobian.parameter_rank == 3
    for observed, expected in zip(result.reporting_rows, expected_reporting, strict=True):
        assert (
            observed.temperature_k,
            observed.predicted_pressure_pa,
            observed.predicted_liquid_density_kg_m3,
        ) == pytest.approx(expected, rel=2.0e-9, abs=2.0e-9)


def test_identity_mismatch_is_rejected_before_native_solve() -> None:
    with pytest.raises(ValueError, match="dataset and specification"):
        fit_pure_saturation(
            model=_model("methane"),
            dataset=load_pure_saturation_dataset("methane"),
            specification=ETHANE_SATURATION_FIT_V1,
        )
    with pytest.raises(ValueError, match="fingerprint"):
        fit_pure_saturation(
            model=_model("methane"),
            dataset=load_pure_saturation_dataset("ethane"),
            specification=ETHANE_SATURATION_FIT_V1,
        )


def test_reporting_conversion_rejects_final_topology_loss() -> None:
    dataset = load_pure_saturation_dataset("methane")
    source = dataset.rows[0]
    native_row = (
        source.row_id,
        source.source_id,
        source.temperature_k,
        source.pressure_pa,
        source.liquid_density_kg_m3,
        source.pressure_pa,
        source.liquid_density_kg_m3,
        9.995e-5,
        1.0e-4,
        1.0,
        1.0,
        (0.0, 0.0, 0.0),
        "CONVERGENCE",
        True,
        "",
    )

    diagnostic = workflow._reporting_row_diagnostic(
        source,
        frozenset(dataset.training_row_ids),
        frozenset(row.row_id for row in dataset.held_out_rows),
        frozenset(),
        METHANE_SATURATION_FIT_V1,
        native_row,
    )

    assert not diagnostic.physically_valid
    assert any("topology separation" in reason for reason in diagnostic.failure_reasons)


def test_provider_callback_failure_is_returned_as_structured_fit_evidence(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    capsule, table, callback = _failing_provider_capsule()
    model = SimpleNamespace(parameter_fingerprint=_model("methane").parameter_fingerprint)
    monkeypatch.setattr(workflow, "native_sdk", lambda _model: capsule)

    result = fit_pure_saturation(
        model=model,
        dataset=load_pure_saturation_dataset("methane"),
        specification=METHANE_SATURATION_FIT_V1,
    )

    assert table.model_context == 1
    assert callback
    assert not result.solver_converged
    assert any("synthetic provider domain failure" in reason for reason in result.failure_reasons)
