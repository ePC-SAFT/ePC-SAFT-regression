from __future__ import annotations

from dataclasses import replace
import math

import pytest

import epcsaft_regression
from epcsaft_regression.records import (
    ETHANE_DATA_SHA256,
    ETHANE_PACKAGED_DATA_SHA256,
    ETHANE_SATURATION_FIT_V1,
    METHANE_DATA_SHA256,
    METHANE_PACKAGED_DATA_SHA256,
    METHANE_SATURATION_FIT_V1,
    SourceIdentity,
    SaturationObservation,
    load_pure_saturation_dataset,
)


METHANE_ROWS = (
    (100.0, 34_375.892, 438.88524),
    (110.0, 88_130.038, 424.77725),
    (120.0, 191_430.08, 409.90234),
    (130.0, 367_319.94, 394.03734),
    (140.0, 641_181.43, 376.86505),
    (150.0, 1_039_961.3, 357.89846),
    (160.0, 1_592_078.0, 336.31495),
    (170.0, 2_328_348.8, 310.50203),
    (180.0, 3_285_180.7, 276.22850),
)
ETHANE_ROWS = (
    (100.0, 11.080787, 640.94852),
    (120.0, 352.30167, 618.94997),
    (140.0, 3_813.5564, 596.58156),
    (160.0, 21_405.224, 573.55144),
    (180.0, 78_638.137, 549.50874),
    (200.0, 217_232.94, 523.97698),
    (220.0, 492_046.39, 496.27145),
    (240.0, 966_787.79, 465.30887),
    (260.0, 1_711_835.4, 429.07617),
    (280.0, 2_806_735.8, 382.72712),
)


def test_public_surface_is_pure_saturation_only() -> None:
    assert hasattr(epcsaft_regression, "load_pure_saturation_dataset")
    assert hasattr(epcsaft_regression, "fit_pure_saturation")
    assert hasattr(epcsaft_regression, "PureSaturationFitResult")
    for retired in (
        "load_methane_dataset",
        "fit_methane_saturation",
        "MethaneFitResult",
        "METHANE_FIT_SPECIFICATION_V1",
    ):
        assert not hasattr(epcsaft_regression, retired)


@pytest.mark.parametrize(
    ("component_id", "expected_rows", "training", "held_out", "stress", "raw_hash", "packaged_hash"),
    (
        (
            "methane",
            METHANE_ROWS,
            (110.0, 130.0, 150.0, 170.0),
            (100.0, 120.0, 140.0, 160.0, 180.0),
            (),
            METHANE_DATA_SHA256,
            METHANE_PACKAGED_DATA_SHA256,
        ),
        (
            "ethane",
            ETHANE_ROWS,
            (140.0, 180.0, 220.0, 260.0),
            (120.0, 160.0, 200.0, 240.0),
            (100.0, 280.0),
            ETHANE_DATA_SHA256,
            ETHANE_PACKAGED_DATA_SHA256,
        ),
    ),
)
def test_retained_dataset_has_exact_source_rows_and_partition(
    component_id: str,
    expected_rows: tuple[tuple[float, float, float], ...],
    training: tuple[float, ...],
    held_out: tuple[float, ...],
    stress: tuple[float, ...],
    raw_hash: str,
    packaged_hash: str,
) -> None:
    dataset = load_pure_saturation_dataset(component_id)

    assert dataset.component_id == component_id
    assert dataset.temperature_unit == "K"
    assert dataset.pressure_unit == "Pa"
    assert dataset.liquid_density_unit == "kg/m3"
    assert dataset.source.retrieved_on == "2026-07-17"
    assert "CRLF line endings normalized to LF" in dataset.source.transformation
    assert dataset.source.data_sha256 == raw_hash
    assert dataset.source.packaged_data_sha256 == packaged_hash
    assert tuple(
        (row.temperature_k, row.pressure_pa, row.liquid_density_kg_m3)
        for row in dataset.rows
    ) == expected_rows
    assert dataset.training_temperatures_k == training
    assert dataset.held_out_temperatures_k == held_out
    assert dataset.stress_temperatures_k == stress


@pytest.mark.parametrize(
    ("component_id", "specification", "start", "molar_mass", "fingerprint"),
    (
        (
            "methane",
            METHANE_SATURATION_FIT_V1,
            (1.08, 3.555744, 157.5315),
            0.016043,
            "sha256:5f836aa84935df70be2e5cffae51b178a7b797c2cee036e9ff47d8097ca94bbf",
        ),
        (
            "ethane",
            ETHANE_SATURATION_FIT_V1,
            (1.6069, 3.5206, 191.42),
            0.030070,
            "sha256:288fbcaa1304881c16f64c3a784eeed19b75c58cca4558f92a21268e5e91258a",
        ),
    ),
)
def test_component_specification_is_explicit_and_dimensionally_fixed(
    component_id: str,
    specification: object,
    start: tuple[float, float, float],
    molar_mass: float,
    fingerprint: str,
) -> None:
    assert specification.component_id == component_id
    assert specification.parameter_names == (
        "segment_count",
        "segment_diameter_angstrom",
        "dispersion_energy_over_k_kelvin",
    )
    assert specification.parameter_units == ("1", "angstrom", "K")
    assert specification.start == start
    assert specification.lower_bounds == (0.5, 2.0, 50.0)
    assert specification.upper_bounds == (3.5, 5.0, 400.0)
    assert specification.parameter_scales == (0.1, 0.1, 10.0)
    assert specification.fixed_amount_mol == 1.0
    assert specification.molar_mass_kg_per_mol == molar_mass
    assert specification.expected_provider_fingerprint == fingerprint
    assert specification.residual_weights == (0.25, 0.25, 0.25, 0.25)
    assert specification.confirmation_parameter_scaled_max_delta == 1.0e-5
    assert specification.confirmation_cost_relative_delta == 1.0e-8
    assert specification.reporting_pressure_scaled_residual_max == 1.0e-8
    assert specification.reporting_chemical_potential_residual_max == 1.0e-8


@pytest.mark.parametrize("value", ("Methane", "ETHANE", "propane", ""))
def test_loader_rejects_aliases_case_variants_and_unknown_strings(value: str) -> None:
    with pytest.raises(ValueError, match="'methane' or 'ethane'"):
        load_pure_saturation_dataset(value)


@pytest.mark.parametrize("value", (None, 1, True, b"methane"))
def test_loader_rejects_non_strings(value: object) -> None:
    with pytest.raises(TypeError, match="exact string"):
        load_pure_saturation_dataset(value)  # type: ignore[arg-type]


@pytest.mark.parametrize(
    ("field", "value", "match"),
    (
        ("row_id", "", "row_id"),
        ("component_id", "propane", "component_id"),
        ("temperature_k", math.nan, "finite"),
        ("pressure_pa", 0.0, "positive"),
        ("liquid_density_kg_m3", math.inf, "finite"),
        ("source_id", "", "source_id"),
    ),
)
def test_observation_rejects_invalid_or_unsupported_fields(
    field: str, value: object, match: str
) -> None:
    valid = SaturationObservation(
        row_id="nist-methane-sat-110-k",
        component_id="methane",
        temperature_k=110.0,
        pressure_pa=88_130.038,
        liquid_density_kg_m3=424.77725,
        source_id="nist-webbook-srd69-methane-saturation",
    )

    with pytest.raises(ValueError, match=match):
        replace(valid, **{field: value})


def test_source_and_dataset_reject_identity_or_partition_changes() -> None:
    dataset = load_pure_saturation_dataset("ethane")

    with pytest.raises(ValueError, match="source fields"):
        replace(dataset.source, packaged_data_sha256="0" * 64)
    with pytest.raises(ValueError, match="source fields"):
        replace(dataset.source, url="https://example.invalid/not-the-admitted-source")
    with pytest.raises(ValueError, match="duplicate row_id"):
        replace(dataset, rows=dataset.rows + (dataset.rows[-1],))
    with pytest.raises(ValueError, match="strictly increasing"):
        replace(dataset, rows=(dataset.rows[1], dataset.rows[0], *dataset.rows[2:]))
    with pytest.raises(ValueError, match="partition"):
        replace(dataset, held_out_temperatures_k=dataset.held_out_temperatures_k[:-1])


def test_source_record_is_frozen() -> None:
    source: SourceIdentity = load_pure_saturation_dataset("methane").source

    with pytest.raises(AttributeError):
        source.source_id = "changed"  # type: ignore[misc]
