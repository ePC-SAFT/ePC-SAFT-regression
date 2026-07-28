from __future__ import annotations

from dataclasses import replace
import math

import pytest

from epcsaft_regression.parameter_regression import (
    AqueousKijMeanIonicActivityObservation,
    AffineParameterTransform,
    ComponentParameterIdentity,
    FixedCompositionVleObservation,
    MeanIonicActivityObservation,
    ModelParameterIdentity,
    ObservationPartition,
    PairParameterIdentity,
    ParameterCoordinate,
    ParameterFamily,
    RegressionProblem,
    RelativePermittivityObservation,
    SolvationGibbsObservation,
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
        (ParameterFamily.BORN_DIAMETER, "angstrom"),
        (ParameterFamily.SOLVATION_FACTOR, "1"),
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


def test_direct_observations_bind_model_order_active_component_and_units() -> None:
    miac = MeanIonicActivityObservation(
        row_id="nabr-001",
        source_id="doi:example",
        source_locator="table-1:nabr-001",
        component_ids=("water", "sodium-cation", "bromide-anion"),
        active_component_id="water",
        temperature_k=298.15,
        pressure_pa=100_000.0,
        formula_unit_molality_mol_per_kg=0.1,
        observed_mean_ionic_activity_coefficient=0.778,
        relative_residual_scale=1.0,
        partition=ObservationPartition.TRAINING,
    )
    gibbs = SolvationGibbsObservation(
        row_id="sodium-s5",
        source_id="doi:example",
        source_locator="table-s5:sodium",
        component_ids=("water", "sodium-cation", "chloride-anion"),
        active_component_id="sodium-cation",
        temperature_k=298.15,
        pressure_pa=100_000.0,
        observed_solvation_gibbs_j_per_mol=-381_100.0,
        residual_scale_j_per_mol=381_100.0,
        partition=ObservationPartition.TRAINING,
    )

    assert miac.component_ids[0] == miac.active_component_id
    assert gibbs.component_ids[1] == gibbs.active_component_id
    assert canonical_dataset_sha256((miac, gibbs))


def test_aqueous_kij_observation_binds_active_pair_and_fixed_context() -> None:
    row = AqueousKijMeanIonicActivityObservation(
        row_id="nabr-001",
        source_id="doi:example",
        source_locator="table-1:nabr-001",
        component_ids=("water", "sodium-cation", "bromide-anion"),
        active_pair_component_ids=("sodium-cation", "water"),
        fixed_k_ij=(-0.3, -0.3, 0.65),
        temperature_k=298.15,
        pressure_pa=100_000.0,
        formula_unit_molality_mol_per_kg=0.1,
        observed_mean_ionic_activity_coefficient=0.778,
        relative_residual_scale=1.0,
        partition=ObservationPartition.TRAINING,
    )

    assert row.active_pair_component_ids == ("sodium-cation", "water")
    assert row.canonical_active_pair_component_ids == ("sodium-cation", "water")
    assert canonical_dataset_sha256((row,))


def test_relative_permittivity_observation_binds_one_model_parameter() -> None:
    row = RelativePermittivityObservation(
        row_id="figiel-water-001",
        source_id="doi:example",
        source_locator="figure-s1:water:001",
        component_ids=("water", "sodium-cation", "bromide-anion"),
        temperature_k=298.15,
        pressure_pa=100_000.0,
        total_ion_mole_fraction=0.05,
        observed_relative_permittivity=58.4,
        residual_scale=1.0,
        partition=ObservationPartition.TRAINING,
    )
    coordinate = ParameterCoordinate(
        family=ParameterFamily.DIELECTRIC_ION_SUPPRESSION_COEFFICIENT,
        identity=ModelParameterIdentity(),
        capability_id="figiel_dielectric_suppression_v1",
        provider_parameter_fingerprint=PARAMETER_HASH,
        provider_topology_fingerprint=TOPOLOGY_HASH,
        unit="1",
        transform=AffineParameterTransform(origin=7.0, scale=1.0),
        lower_bound=0.01,
        upper_bound=30.0,
        starts=(2.0, 12.0),
    )

    assert coordinate.identity.canonical_component_ids == ()
    assert row.component_ids == (
        "water",
        "sodium-cation",
        "bromide-anion",
    )
    assert canonical_dataset_sha256((row,))


@pytest.mark.parametrize(
    ("build", "match"),
    (
        (
            lambda: MeanIonicActivityObservation(
                row_id="nabr-001",
                source_id="doi:example",
                source_locator="table-1:nabr-001",
                component_ids=("water", "sodium-cation", "bromide-anion"),
                active_component_id="methane",
                temperature_k=298.15,
                pressure_pa=100_000.0,
                formula_unit_molality_mol_per_kg=0.1,
                observed_mean_ionic_activity_coefficient=0.778,
                relative_residual_scale=1.0,
                partition=ObservationPartition.TRAINING,
            ),
            "active_component_id",
        ),
        (
            lambda: SolvationGibbsObservation(
                row_id="sodium-s5",
                source_id="doi:example",
                source_locator="table-s5:sodium",
                component_ids=("water", "sodium-cation", "chloride-anion"),
                active_component_id="sodium-cation",
                temperature_k=298.15,
                pressure_pa=100_000.0,
                observed_solvation_gibbs_j_per_mol=-381_100.0,
                residual_scale_j_per_mol=0.0,
                partition=ObservationPartition.TRAINING,
            ),
            "positive",
        ),
    ),
)
def test_direct_observations_reject_invalid_identity_or_scale(
    build: object, match: str
) -> None:
    with pytest.raises((TypeError, ValueError), match=match):
        build()  # type: ignore[operator]


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
