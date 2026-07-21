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
    PROPANE_FIT_TARGET_CONTRACT_SHA256,
    PROPANE_PACKAGED_DATA_SHA256,
    PROPANE_PACKET_YAML_SHA256,
    PROPANE_SATURATION_FIT_V1,
    PROPANE_SOURCE_RECEIPT_SHA256,
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
PROPANE_ROWS = (
    (110.0, 0.6, 707.968),
    (120.0, 3.2, 697.825),
    (130.0, 18.0, 687.713),
    (140.0, 78.0, 677.601),
    (150.0, 283.0, 667.462),
    (160.0, 851.0, 657.272),
    (170.0, 2_205.0, 647.004),
    (180.0, 5_068.0, 636.628),
    (190.0, 10_547.0, 626.123),
    (200.0, 20_193.0, 615.456),
    (210.0, 36_032.0, 604.594),
    (220.0, 60_574.0, 593.499),
    (230.0, 96_775.0, 582.132),
    (240.0, 148_000.0, 570.444),
    (250.0, 217_964.0, 558.383),
    (260.0, 310_685.0, 545.88),
    (270.0, 430_425.0, 532.853),
    (280.0, 581_684.0, 519.198),
    (290.0, 769_143.0, 504.796),
    (300.0, 997_682.0, 489.465),
    (310.0, 1_272_430.0, 472.968),
    (320.0, 1_598_870.0, 454.951),
    (330.0, 1_983_000.0, 434.869),
    (340.0, 2_431_450.0, 411.772),
)


def test_public_surface_preserves_pure_saturation_without_retired_aliases() -> None:
    assert hasattr(epcsaft_regression, "load_pure_saturation_dataset")
    assert hasattr(epcsaft_regression, "fit_pure_saturation")
    assert hasattr(epcsaft_regression, "PureSaturationFitResult")
    assert epcsaft_regression.PROPANE_SATURATION_FIT_V1 is PROPANE_SATURATION_FIT_V1
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
        (
            "propane",
            PROPANE_ROWS,
            (150.0, 210.0, 270.0, 330.0),
            tuple(float(value) for value in range(120, 330, 10) if value not in (150, 210, 270)),
            (110.0, 340.0),
            PROPANE_PACKAGED_DATA_SHA256,
            PROPANE_PACKAGED_DATA_SHA256,
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
    if component_id in ("methane", "ethane"):
        assert "CRLF line endings normalized to LF" in dataset.source.transformation
    else:
        assert "Validation commit 7e51590757f1cb85f51df98e9fe1f88cd4255a88" in (
            dataset.source.transformation
        )
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
        (
            "propane",
            PROPANE_SATURATION_FIT_V1,
            (2.002, 3.6184, 208.11),
            0.044096,
            "sha256:9bfbc8d7789e51609945e61dbdf7a020decc8f9e31b408b0977724c7cb3e1551",
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
    assert specification.max_num_iterations == (5000 if component_id == "propane" else 500)


@pytest.mark.parametrize("value", ("Methane", "ETHANE", "Propane", "PROPANE", "butane", ""))
def test_loader_rejects_aliases_case_variants_and_unknown_strings(value: str) -> None:
    with pytest.raises(ValueError, match="'methane', 'ethane', or 'propane'"):
        load_pure_saturation_dataset(value)


@pytest.mark.parametrize("value", (None, 1, True, b"methane"))
def test_loader_rejects_non_strings(value: object) -> None:
    with pytest.raises(TypeError, match="exact string"):
        load_pure_saturation_dataset(value)  # type: ignore[arg-type]


@pytest.mark.parametrize(
    ("field", "value", "match"),
    (
        ("row_id", "", "row_id"),
        ("component_id", "butane", "component_id"),
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


def test_propane_packet_identity_and_uncertainties_are_retained_without_cutoffs() -> None:
    dataset = load_pure_saturation_dataset("propane")

    assert PROPANE_PACKAGED_DATA_SHA256 == (
        "ccd1cfa15ec44432b06cbf22316d168c61b282631c9b1e1591e497b8d48b5676"
    )
    assert PROPANE_PACKET_YAML_SHA256 == (
        "ba31448989f565d05d63908076e836977780aa87199f208310e9b80b03f64697"
    )
    assert PROPANE_SOURCE_RECEIPT_SHA256 == (
        "ed5eb703ccd3e6bb4c4cfa82ecd58c58f9da0c93ab07a204dee94d8b0ae8d081"
    )
    assert PROPANE_FIT_TARGET_CONTRACT_SHA256 == (
        "7f25259265dfa42f1de36bc04740baf6c78e09c8bc35a42392f06a4b8a32cb90"
    )
    assert tuple(
        (
            row.pressure_expanded_uncertainty_pa,
            row.liquid_density_expanded_uncertainty_kg_m3,
            row.vapor_density_kg_m3,
            row.vapor_density_expanded_uncertainty_kg_m3,
        )
        for row in dataset.rows
    ) == (
        (0.6, 0.133, None, None),
        (1.6, 0.131, None, None),
        (1.1, 0.13, None, None),
        (3.1, 0.128, None, None),
        (5.0, 0.127, None, None),
        (9.0, 0.125, None, None),
        (10.0, 0.124, None, None),
        (11.0, 0.122, None, None),
        (16.0, 0.121, None, None),
        (19.0, 0.119, None, None),
        (23.0, 0.118, None, None),
        (19.0, 0.117, None, None),
        (24.0, 0.116, 2.3148, 0.0051),
        (33.0, 0.115, 3.4379, 0.0063),
        (37.0, 0.114, 4.9398, 0.0081),
        (50.0, 0.113, 6.9029, 0.0098),
        (65.0, 0.113, 9.423, 0.0125),
        (70.0, 0.113, 12.6155, 0.0158),
        (77.0, 0.114, 16.6231, 0.0201),
        (90.0, 0.115, 21.6299, 0.0257),
        (100.0, 0.118, 27.8854, 0.0331),
        (110.0, 0.124, 35.7463, 0.043),
        (140.0, 0.132, 45.7626, 0.057),
        (170.0, 0.138, 58.878, 0.0666),
    )
    assert dataset.source.use_basis.endswith(
        "source evidence, not model-acceptance cutoffs"
    )
