from __future__ import annotations

from dataclasses import replace
import math

import pytest

from epcsaft_regression.parameter_regression import (
    AffineParameterTransform,
    ComponentParameterIdentity,
    FixedCompositionVleObservation,
    ObservationPartition,
    PairParameterIdentity,
    ParameterCoordinate,
    ParameterFamily,
    RegressionProblem,
    SourceDescriptor,
    canonical_dataset_sha256,
)


SOURCE_HASH = "a" * 64
PARAMETER_HASH = "sha256:" + "b" * 64
TOPOLOGY_HASH = "sha256:" + "c" * 64


def _row(
    row_id: str,
    partition: ObservationPartition = ObservationPartition.TRAINING,
) -> FixedCompositionVleObservation:
    return FixedCompositionVleObservation(
        row_id=row_id,
        source_id="doi:example",
        source_locator=f"table-1:{row_id}",
        component_ids=("methane", "ethane"),
        temperature_k=200.0,
        pressure_pa=2.0e6,
        liquid_mole_fraction_first=0.25,
        vapor_mole_fraction_first=0.75,
        pressure_scale_pa=2.0e6,
        chemical_potential_scales=(1.0, 1.0),
        liquid_volume_origin_m3_per_mol=6.0e-5,
        liquid_volume_start_m3_per_mol=6.5e-5,
        liquid_volume_bounds_m3_per_mol=(2.0e-5, 1.0e-4),
        vapor_volume_origin_m3_per_mol=8.0e-4,
        vapor_volume_start_m3_per_mol=9.0e-4,
        vapor_volume_bounds_m3_per_mol=(1.0e-4, 1.0e-2),
        partition=partition,
    )


def _problem(
    rows: tuple[FixedCompositionVleObservation, ...],
    *,
    identity: PairParameterIdentity | None = None,
    dataset_hash: str | None = None,
) -> RegressionProblem:
    source = SourceDescriptor(
        source_id="doi:example",
        citation="Example et al. (2026)",
        durable_locator="https://doi.org/10.0000/example",
        source_artifact_sha256=SOURCE_HASH,
        canonical_dataset_sha256=dataset_hash or canonical_dataset_sha256(rows),
        transformation_record="Converted bar to Pa; no rows omitted.",
        units_and_bases="T/K; P/Pa; liquid and vapor mole fraction.",
        use_basis="Source-backed regression evidence.",
        residual_scale_rationale="Pressure scaled by observed pressure; mu/RT dimensionless.",
    )
    coordinate = ParameterCoordinate(
        family=ParameterFamily.K_IJ,
        identity=identity or PairParameterIdentity("methane", "ethane"),
        capability_id="neutral_binary_phase_kij_v1",
        provider_parameter_fingerprint=PARAMETER_HASH,
        provider_topology_fingerprint=TOPOLOGY_HASH,
        unit="1",
        transform=AffineParameterTransform(origin=0.0, scale=0.01),
        lower_bound=-0.15,
        upper_bound=0.10,
        starts=(0.0, -0.05, 0.05),
    )
    return RegressionProblem(
        sources=(source,),
        parameters=(coordinate,),
        observations=rows,
        maximum_condition_number=1.0e10,
        maximum_iterations=200,
        function_tolerance=1.0e-12,
        gradient_tolerance=1.0e-12,
        parameter_tolerance=1.0e-12,
        confirmation_parameter_scaled_max_delta=1.0e-5,
        confirmation_cost_relative_delta=1.0e-8,
    )


def test_contract_canonicalizes_pair_and_binds_source_dataset_and_partitions() -> None:
    rows = (
        _row("train-1"),
        _row("held-1", ObservationPartition.HELD_OUT),
        _row("stress-1", ObservationPartition.STRESS),
    )
    problem = _problem(rows, identity=PairParameterIdentity("ethane", "methane"))

    assert problem.parameters[0].identity.component_ids == ("ethane", "methane")
    assert problem.parameters[0].identity.canonical_component_ids == ("ethane", "methane")
    assert canonical_dataset_sha256(rows) == problem.sources[0].canonical_dataset_sha256
    assert problem.training_observations == (rows[0],)
    assert problem.held_out_observations == (rows[1],)
    assert problem.stress_observations == (rows[2],)
    assert problem.parameters[0].solver_starts == (0.0, -5.0, 5.0)


def test_pair_coordinate_accepts_the_distinct_lij_family() -> None:
    coordinate = replace(
        _problem((_row("train-1"),)).parameters[0],
        family=ParameterFamily.L_IJ,
        capability_id="neutral_binary_phase_lij_v1",
    )

    assert coordinate.family is ParameterFamily.L_IJ
    assert coordinate.identity.canonical_component_ids == ("ethane", "methane")
    assert coordinate.unit == "1"


@pytest.mark.parametrize(
    ("family", "unit"),
    (
        (ParameterFamily.SEGMENT_COUNT, "1"),
        (ParameterFamily.SEGMENT_DIAMETER, "angstrom"),
        (ParameterFamily.DISPERSION_ENERGY_OVER_K, "K"),
    ),
)
def test_component_coordinate_accepts_scalar_pure_families(
    family: ParameterFamily,
    unit: str,
) -> None:
    base = _problem((_row("train-1"),)).parameters[0]
    coordinate = replace(
        base,
        family=family,
        identity=ComponentParameterIdentity("methane"),
        capability_id=f"neutral_pure_{family.value}_v1",
        unit=unit,
        lower_bound=0.1,
        upper_bound=500.0,
        starts=(1.0, 1.1),
    )

    assert coordinate.identity.canonical_component_ids == ("methane",)
    assert coordinate.unit == unit


@pytest.mark.parametrize(
    ("mutation", "match"),
    (
        (lambda rows: replace(rows[0], row_id=""), "row_id"),
        (lambda rows: replace(rows[0], component_ids=("methane", "methane")), "distinct"),
        (lambda rows: replace(rows[0], temperature_k=math.nan), "finite"),
        (lambda rows: replace(rows[0], pressure_pa=0.0), "positive"),
        (lambda rows: replace(rows[0], liquid_mole_fraction_first=1.0), "strictly between"),
        (lambda rows: replace(rows[0], pressure_scale_pa=math.inf), "finite"),
        (
            lambda rows: replace(rows[0], liquid_volume_start_m3_per_mol=2.0e-4),
            "liquid volume start",
        ),
        (
            lambda rows: replace(rows[0], vapor_volume_bounds_m3_per_mol=(1.0e-2, 1.0e-4)),
            "vapor volume bounds",
        ),
    ),
)
def test_observation_rejects_ambiguous_or_nonphysical_values(mutation: object, match: str) -> None:
    rows = (_row("train-1"),)
    with pytest.raises((TypeError, ValueError), match=match):
        mutation(rows)  # type: ignore[operator]


@pytest.mark.parametrize(
    ("build", "match"),
    (
        (
            lambda rows: _problem(rows, dataset_hash="d" * 64),
            "canonical dataset SHA-256",
        ),
        (
            lambda rows: _problem((rows[0], replace(rows[0], source_locator="other"),)),
            "duplicate row_id",
        ),
        (
            lambda rows: replace(
                _problem(rows).parameters[0],
                identity=PairParameterIdentity("methane", "methane"),
            ),
            "distinct",
        ),
        (
            lambda rows: replace(
                _problem(rows).parameters[0],
                provider_topology_fingerprint="",
            ),
            "topology fingerprint",
        ),
        (
            lambda rows: replace(
                _problem(rows).parameters[0],
                lower_bound=0.10,
                upper_bound=-0.15,
            ),
            "bounds",
        ),
        (
            lambda rows: replace(
                _problem(rows).parameters[0],
                starts=(0.2, 0.3),
            ),
            "start",
        ),
        (
            lambda rows: replace(
                _problem(rows),
                parameters=(
                    _problem(rows).parameters[0],
                    _problem(rows).parameters[0],
                ),
            ),
            "duplicate parameter identity",
        ),
    ),
)
def test_problem_fails_closed_on_identity_hash_bounds_and_partition_errors(
    build: object,
    match: str,
) -> None:
    rows = (_row("train-1"),)
    with pytest.raises((TypeError, ValueError), match=match):
        build(rows)  # type: ignore[operator]
