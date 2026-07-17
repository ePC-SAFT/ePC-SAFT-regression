from __future__ import annotations

from dataclasses import replace
import math

import pytest

from epcsaft_regression.records import (
    EXPECTED_DATA_SHA256,
    METHANE_FIT_SPECIFICATION_V1,
    SourceIdentity,
    SaturationObservation,
    load_methane_dataset,
)


EXPECTED_ROWS = (
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


def test_retained_methane_dataset_has_exact_source_rows_and_training_partition() -> None:
    dataset = load_methane_dataset()

    assert dataset.dataset_id == "nist-webbook-methane-saturation-100-180-k-v1"
    assert dataset.species == "methane"
    assert dataset.temperature_unit == "K"
    assert dataset.pressure_unit == "Pa"
    assert dataset.liquid_density_unit == "kg/m3"
    assert dataset.source.source_id == "nist-webbook-srd69-methane-saturation"
    assert dataset.source.retrieved_on == "2026-07-17"
    assert "CRLF line endings normalized to LF" in dataset.source.transformation
    assert dataset.source.data_sha256 == EXPECTED_DATA_SHA256 == (
        "a5e16df3bf8ec78483fc340782cddc89ab8b284a9f6dfaecd6cda3ffde579227"
    )
    assert tuple(
        (row.temperature_k, row.pressure_pa, row.liquid_density_kg_m3)
        for row in dataset.rows
    ) == EXPECTED_ROWS
    assert dataset.training_row_ids == (
        "nist-methane-sat-110-k",
        "nist-methane-sat-130-k",
        "nist-methane-sat-150-k",
        "nist-methane-sat-170-k",
    )


def test_first_slice_specification_is_explicit_and_dimensionally_fixed() -> None:
    specification = METHANE_FIT_SPECIFICATION_V1

    assert specification.parameter_names == (
        "segment_count",
        "segment_diameter_angstrom",
        "dispersion_energy_over_k_kelvin",
    )
    assert specification.parameter_units == ("1", "angstrom", "K")
    assert specification.start == (1.08, 3.555744, 157.5315)
    assert specification.lower_bounds == (0.5, 2.0, 50.0)
    assert specification.upper_bounds == (3.5, 5.0, 400.0)
    assert specification.parameter_scales == (0.1, 0.1, 10.0)
    assert specification.fixed_amount_mol == 1.0
    assert specification.methane_molar_mass_kg_per_mol == 0.016043
    assert specification.residual_names == (
        "liquid_pressure",
        "vapor_pressure",
        "chemical_potential_equality",
        "liquid_density",
    )
    assert specification.residual_weights == (0.25, 0.25, 0.25, 0.25)
    assert specification.liquid_volume_bounds_m3 == (2.0e-5, 1.0e-4)
    assert specification.vapor_volume_bounds_m3 == (1.5e-4, 0.1)
    assert specification.training_temperatures_k == (110.0, 130.0, 150.0, 170.0)
    assert specification.reporting_pressure_bounds_pa == (1.0e3, 1.0e7)
    assert specification.confirmation_liquid_volume_start_multiplier == 1.01
    assert specification.confirmation_vapor_volume_start_multiplier == 0.98
    assert specification.confirmation_parameter_scaled_max_delta == 1.0e-5
    assert specification.confirmation_cost_relative_delta == 1.0e-8
    assert specification.reporting_pressure_scaled_residual_max == 1.0e-8
    assert specification.reporting_chemical_potential_residual_max == 1.0e-8
    assert specification.ceres_linear_solver == "DENSE_QR"
    assert specification.ceres_num_threads == 1
    assert specification.ceres_logging == "SILENT"


@pytest.mark.parametrize(
    ("field", "value", "match"),
    (
        ("row_id", "", "row_id"),
        ("species", "ethane", "methane"),
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
        species="methane",
        temperature_k=110.0,
        pressure_pa=88_130.038,
        liquid_density_kg_m3=424.77725,
        source_id="nist-webbook-srd69-methane-saturation",
    )

    with pytest.raises(ValueError, match=match):
        replace(valid, **{field: value})


def test_source_identity_rejects_missing_units_or_changed_hash() -> None:
    dataset = load_methane_dataset()

    with pytest.raises(ValueError, match="units"):
        replace(dataset.source, units=())
    with pytest.raises(ValueError, match="SHA-256"):
        replace(dataset.source, data_sha256="0" * 64)
    with pytest.raises(ValueError, match="exact retained source identity"):
        replace(dataset.source, url="https://example.invalid/not-the-admitted-source")
    with pytest.raises(ValueError, match="exact retained source identity"):
        replace(dataset.source, citation="arbitrary citation")


def test_dataset_rejects_duplicates_wrong_order_and_wrong_training_rows() -> None:
    dataset = load_methane_dataset()

    with pytest.raises(ValueError, match="duplicate row_id"):
        replace(dataset, rows=dataset.rows + (dataset.rows[-1],))
    with pytest.raises(ValueError, match="strictly increasing"):
        replace(dataset, rows=(dataset.rows[1], dataset.rows[0], *dataset.rows[2:]))
    with pytest.raises(ValueError, match="training partition"):
        replace(dataset, training_row_ids=dataset.training_row_ids[:-1])
    swapped_ids = (
        replace(dataset.rows[0], row_id=dataset.rows[1].row_id),
        replace(dataset.rows[1], row_id=dataset.rows[0].row_id),
        *dataset.rows[2:],
    )
    with pytest.raises(ValueError, match="row IDs must match"):
        replace(dataset, rows=swapped_ids)


def test_source_record_is_frozen() -> None:
    source: SourceIdentity = load_methane_dataset().source

    with pytest.raises(AttributeError):
        source.source_id = "changed"  # type: ignore[misc]
