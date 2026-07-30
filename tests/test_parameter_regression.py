from __future__ import annotations

import csv
import math
from dataclasses import asdict, replace
from pathlib import Path

from epcsaft import Mixture, Parameters, native_sdk, unit_registry
from epcsaft.records import (
    AssociationParameterRecord,
    ComponentRecord,
    ModelParameterRecord,
    PairParameterRecord,
    SingleParameterRecord,
    SiteRecord,
    SourceRecord,
    ValidityDomain,
)
import pytest

import epcsaft_regression.parameter_regression as parameter_regression
from epcsaft_regression import (
    AcquisitionClass,
    AffineParameterTransform,
    AqueousKijMeanIonicActivityObservation,
    ComponentParameterIdentity,
    ConfirmationControls,
    DirectObservationRowDiagnostic,
    FixedCompositionVleObservation,
    FIGIEL_BORN_DIAMETER_TRACER_V1,
    FIGIEL_AQUEOUS_KIJ_V1,
    FIGIEL_WATER_SOLVATION_FACTOR_V1,
    IonSolvationKijObservation,
    MeanIonicActivityObservation,
    ModelParameterIdentity,
    ObservationPartition,
    ObjectiveContract,
    ObservationDataset,
    PairParameterIdentity,
    ParameterCoordinate,
    ParameterFamily,
    ParameterRequest,
    PureDensityObservation,
    PureDensityRowDiagnostic,
    PureSaturationObservation,
    PureSaturationRowDiagnostic,
    RegressionProblem,
    RegressionResult,
    RankControls,
    RelativePermittivityRatioObservation,
    SourceDescriptor,
    SourceInput,
    RowProvenance,
    SolverControls,
    SolvationGibbsObservation,
    UnsupportedParameterCapability,
    canonical_dataset_sha256,
    fit_parameters,
    load_pure_saturation_dataset,
    parameter_capabilities,
    prepare_fit,
)
from epcsaft_regression.workflow import _aqueous_kij_models, _fixed_water_factor_model
from epcsaft_regression.parameter_regression import (
    _evaluate_parameters,
    _native_payload,
)

_MAY_METHANE_PROPANE_CANONICAL_SHA256 = (
    "c7506ce654d9b6df60ec7ff6bdc6dfde526f82a3d69fc524ac9186976785cefe"
)


def _model() -> Mixture:
    parameters = Parameters.from_catalog(
        "gross-2001-methane-ethane",
        components=("methane", "ethane"),
        version=1,
    )
    return Mixture(parameters)


def _methane_propane_parameters() -> Parameters:
    """Build the campaign's public user-provided two-component parameters."""
    methane_catalog = Parameters.from_catalog(
        "gross-2001-methane-ethane",
        components=("methane",),
        version=1,
    )
    propane_catalog = Parameters.from_catalog(
        "gross-2001-propane",
        components=("propane",),
        version=1,
    )
    methane_records = tuple(
        record
        for record in methane_catalog.records
        if isinstance(record, SingleParameterRecord)
        and record.component_id == "methane"
    )
    propane_records = tuple(
        record
        for record in propane_catalog.records
        if isinstance(record, SingleParameterRecord)
        and record.component_id == "propane"
    )
    source = SourceRecord(
        "gross-sadowski-2001",
        (
            "Gross, J.; Sadowski, G. PC-SAFT: An Equation of State Based on a "
            "Perturbation Theory for Chain Molecules. Industrial & Engineering "
            "Chemistry Research 2001, 40, 1244-1260."
        ),
        "source-backed installed parameter catalog",
        "10.1021/ie0003887",
    )
    domains = tuple(
        {
            domain.domain_id: domain
            for domain in (*methane_catalog.domains, *propane_catalog.domains)
        }.values()
    )
    components = tuple(
        component
        for component in (*methane_catalog.components, *propane_catalog.components)
        if component.component_id in {"methane", "propane"}
    )
    pair = PairParameterRecord(
        "methane-propane-kij-initial",
        "methane",
        "propane",
        "k_ij",
        0.0,
        source.source_id,
        "test-only zero active-pair initialization",
        "gross-pair-unknown",
    )
    return Parameters.from_records(
        bundle_id="may-2015-methane-propane-kij-fixture",
        bundle_version=1,
        purpose="user-provided",
        sources=(source,),
        domains=domains,
        components=components,
        singles=(*methane_records, *propane_records),
        pairs=(pair,),
        models=(
            record
            for record in methane_catalog.records
            if isinstance(record, ModelParameterRecord)
        ),
        selected_components=("methane", "propane"),
    )


def _methane_propane_model() -> Mixture:
    return Mixture(_methane_propane_parameters())


def _pure_model() -> Mixture:
    parameters = Parameters.from_catalog(
        "gross-2001-methane-ethane",
        components=("methane",),
        version=1,
    )
    return Mixture(parameters)


def _associating_pure_model() -> Mixture:
    parameters = Parameters.from_catalog(
        "figiel-2025-reference-electrolytes",
        components=("ethanol",),
        version=1,
    )
    return Mixture(parameters)


def _generic_associating_model(
    sites: tuple[tuple[str, int], ...],
    pairs: tuple[tuple[str, str, float, float], ...],
    *,
    segment_count: float = 3.2,
    segment_diameter_angstrom: float = 3.5,
    dispersion_energy_over_k: float = 280.0,
) -> Mixture:
    provenance = {
        "source_id": "manufactured-association",
        "locator": "generic Regression contract test",
        "domain_id": "manufactured-association-domain",
    }
    associations = tuple(
        record
        for left, right, energy, volume in pairs
        for record in (
            AssociationParameterRecord(
                f"{left}-{right}-energy",
                "test-amine",
                left,
                "test-amine",
                right,
                "association_energy_over_k",
                energy * unit_registry.kelvin,
                **provenance,
            ),
            AssociationParameterRecord(
                f"{left}-{right}-volume",
                "test-amine",
                left,
                "test-amine",
                right,
                "association_volume",
                volume,
                **provenance,
            ),
        )
    )
    return Mixture(
        Parameters.from_records(
            bundle_id="manufactured-generic-association",
            bundle_version=1,
            purpose="user-provided",
            sources=(
                SourceRecord(
                    "manufactured-association",
                    "Manufactured generic association model",
                    "Regression contract test",
                ),
            ),
            domains=(
                ValidityDomain(
                    "manufactured-association-domain",
                    "reported-conditions",
                    temperature_min=250.0 * unit_registry.kelvin,
                    temperature_max=450.0 * unit_registry.kelvin,
                    pressure_min=1.0 * unit_registry.pascal,
                    pressure_max=10.0 * unit_registry.megapascal,
                ),
            ),
            components=(ComponentRecord("test-amine"),),
            singles=(
                SingleParameterRecord(
                    "test-m",
                    "test-amine",
                    "segment_count",
                    segment_count,
                    **provenance,
                ),
                SingleParameterRecord(
                    "test-sigma",
                    "test-amine",
                    "segment_diameter",
                    segment_diameter_angstrom * unit_registry.angstrom,
                    **provenance,
                ),
                SingleParameterRecord(
                    "test-epsilon",
                    "test-amine",
                    "dispersion_energy_over_k",
                    dispersion_energy_over_k * unit_registry.kelvin,
                    **provenance,
                ),
                SingleParameterRecord(
                    "test-molar-mass",
                    "test-amine",
                    "molar_mass",
                    0.088 * unit_registry.kilogram / unit_registry.mole,
                    **provenance,
                ),
            ),
            sites=tuple(
                SiteRecord(
                    f"test-site-{site}",
                    "test-amine",
                    site,
                    site,
                    multiplicity,
                    **provenance,
                )
                for site, multiplicity in sites
            ),
            associations=associations,
            models=(
                ModelParameterRecord(
                    "test-permittivity",
                    "relative_permittivity_formulation",
                    "none",
                    **provenance,
                ),
            ),
            selected_components=("test-amine",),
        )
    )


def _generic_associating_problem(
    model: Mixture,
    pairs: tuple[tuple[str, str, float, float], ...],
) -> RegressionProblem:
    capability = next(
        capability
        for capability in parameter_capabilities(model)
        if not isinstance(capability, UnsupportedParameterCapability)
        and capability.capability_id
        == "neutral_pure_associating_joint_sigma_basis_v1"
    )
    physical = (
        3.2,
        3.5,
        280.0,
        *(value for pair in pairs for value in pair[2:]),
    )
    families = (
        ParameterFamily.SEGMENT_COUNT,
        ParameterFamily.SEGMENT_DIAMETER,
        ParameterFamily.DISPERSION_ENERGY_OVER_K,
        *(
            family
            for _ in pairs
            for family in (
                ParameterFamily.ASSOCIATION_ENERGY_OVER_K,
                ParameterFamily.ASSOCIATION_VOLUME,
            )
        ),
    )
    identities = (
        ComponentParameterIdentity("test-amine"),
        ComponentParameterIdentity("test-amine"),
        ComponentParameterIdentity("test-amine"),
        ModelParameterIdentity(),
        ModelParameterIdentity(),
    )
    units = (
        "1",
        "angstrom",
        "K",
        *(unit for _ in pairs for unit in ("K", "1")),
    )
    scales = (
        0.5,
        0.2,
        50.0,
        *(scale for _ in pairs for scale in (500.0, 0.05)),
    )
    bounds = (
        (1.0, 8.0),
        (2.0, 6.0),
        (50.0, 800.0),
        *(bound for _ in pairs for bound in ((100.0, 8000.0), (1e-5, 0.2))),
    )
    row = PureDensityObservation(
        row_id="manufactured-density",
        source_id="manufactured-observation",
        source_locator="manufactured generic association row",
        component_id="test-amine",
        temperature_k=350.0,
        pressure_pa=100_000.0,
        density_kg_per_m3=800.0,
        molar_mass_kg_per_mol=0.088,
        pressure_scale_pa=100_000.0,
        density_scale_kg_per_m3=800.0,
        volume_origin_m3_per_mol=0.088 / 800.0,
        volume_start_m3_per_mol=0.088 / 800.0,
        volume_bounds_m3_per_mol=(5e-5, 2e-4),
        partition=ObservationPartition.TRAINING,
    )
    source = SourceDescriptor(
        source_id=row.source_id,
        citation="Manufactured generic association row",
        durable_locator=row.source_locator,
        source_artifact_sha256="0" * 64,
        canonical_dataset_sha256=canonical_dataset_sha256((row,)),
        transformation_record="none",
        units_and_bases="SI",
        use_basis="fixed ordinary-sigma neutral-pure-2B regression contract",
        residual_scale_rationale="manufactured finite scales",
    )
    return RegressionProblem(
        sources=(source,),
        parameters=tuple(
            ParameterCoordinate(
                family=family,
                identity=identity,
                capability_id=capability.capability_id,
                provider_parameter_fingerprint=capability.parameter_fingerprint,
                provider_topology_fingerprint=capability.topology_fingerprint,
                unit=unit,
                transform=AffineParameterTransform(origin=value, scale=scale),
                lower_bound=bound[0],
                upper_bound=bound[1],
            )
            for family, identity, unit, value, scale, bound in zip(
                families, identities, units, physical, scales, bounds, strict=True
            )
        ),
        parameter_slot_indices=tuple(range(len(physical))),
        start_vectors=(tuple(physical), tuple(physical)),
        observations=(row,),
        maximum_condition_number=1e12,
        maximum_iterations=20,
        maximum_solver_time_seconds=5.0,
        function_tolerance=1e-10,
        gradient_tolerance=1e-10,
        parameter_tolerance=1e-10,
        confirmation_parameter_scaled_max_delta=1e-6,
        confirmation_cost_relative_delta=1e-6,
    )


def _pure_density_problem(
    model: Mixture,
    family: ParameterFamily,
) -> RegressionProblem:
    capability = _capability(model, family)
    row = PureDensityObservation(
        row_id="held-2012-ethanol-density-298-15-k",
        source_id="held-2012-pure-ethanol-density",
        source_locator=(
            "Validation data/held-2012-pure-ethanol-density.csv:2; "
            "DOI 10.1016/j.ces.2011.09.040, section 2.2"
        ),
        component_id="ethanol",
        temperature_k=298.15,
        pressure_pa=100_000.0,
        density_kg_per_m3=785.54,
        molar_mass_kg_per_mol=0.046069,
        pressure_scale_pa=100_000.0,
        density_scale_kg_per_m3=785.54,
        volume_origin_m3_per_mol=0.046069 / 785.54,
        volume_start_m3_per_mol=0.046069 / 785.54,
        volume_bounds_m3_per_mol=(4.0e-5, 1.0e-4),
        partition=ObservationPartition.TRAINING,
    )
    if family is ParameterFamily.ASSOCIATION_ENERGY_OVER_K:
        origin, scale, bounds, starts, unit = (
            2653.4,
            500.0,
            (1000.0, 5000.0),
            (2200.0, 3200.0),
            "K",
        )
    else:
        origin, scale, bounds, starts, unit = (
            0.03238,
            0.01,
            (0.005, 0.1),
            (0.02, 0.05),
            "1",
        )
    source = SourceDescriptor(
        source_id=row.source_id,
        citation=(
            "Held, Prinz, Wallmeyer, and Sadowski (2012), Chemical "
            "Engineering Science 68, 328-339, DOI 10.1016/j.ces.2011.09.040"
        ),
        durable_locator=row.source_locator,
        source_artifact_sha256=(
            "25e3be94ee3cfb5eb13df89827f4368673f87f96d6b0225c12e35f9396b8779c"
        ),
        canonical_dataset_sha256=canonical_dataset_sha256((row,)),
        transformation_record=(
            "The SI density is transcribed unchanged. The source reports "
            "ambient pressure; the model evaluation uses the predeclared "
            "Validation convention of 1 bar."
        ),
        units_and_bases=(
            "T/K and mass density/(kg/m3) from the source; model pressure/Pa."
        ),
        use_basis=(
            "Direct experimental pure-liquid density anchor with 1 bar as "
            "the explicit model-pressure approximation."
        ),
        residual_scale_rationale=(
            "Observed pressure and density magnitudes are numerical scales, "
            "not experimental uncertainties."
        ),
    )
    return RegressionProblem(
        sources=(source,),
        parameters=(
            ParameterCoordinate(
                family=family,
                identity=ModelParameterIdentity(),
                capability_id=capability.capability_id,
                provider_parameter_fingerprint=capability.parameter_fingerprint,
                provider_topology_fingerprint=capability.topology_fingerprint,
                unit=unit,
                transform=AffineParameterTransform(origin=origin, scale=scale),
                lower_bound=bounds[0],
                upper_bound=bounds[1],
            ),
        ),
        parameter_slot_indices=(0,),
        start_vectors=tuple((start,) for start in starts),
        observations=(row,),
        maximum_condition_number=1.0e12,
        maximum_iterations=100,
        maximum_solver_time_seconds=30.0,
        function_tolerance=1.0e-12,
        gradient_tolerance=1.0e-12,
        parameter_tolerance=1.0e-12,
        confirmation_parameter_scaled_max_delta=1.0e-7,
        confirmation_cost_relative_delta=1.0e-8,
    )


def _aqueous_model(
    component_ids: tuple[str, str, str] = (
        "water",
        "sodium-cation",
        "bromide-anion",
    ),
) -> Mixture:
    parameters = Parameters.from_catalog(
        "figiel-2025-reference-electrolytes",
        components=component_ids,
        version=1,
    )
    return Mixture(parameters)


def _capability(
    model: Mixture,
    family: ParameterFamily = ParameterFamily.K_IJ,
):
    return next(
        capability
        for capability in parameter_capabilities(model)
        if not isinstance(capability, UnsupportedParameterCapability)
        and capability.family is family
    )


def _prepare_existing_problem(
    model: Mixture,
    problem: RegressionProblem,
    residual_family: str,
):
    observation_type = type(problem.observations[0])
    assert all(type(row) is observation_type for row in problem.observations)
    source = problem.sources[0]
    dataset = ObservationDataset.from_records(
        observation_type,
        tuple(asdict(row) for row in problem.observations),
        source=SourceInput(
            source.source_id,
            source.citation,
            source.durable_locator,
            source.source_artifact_sha256,
            source.transformation_record,
            source.units_and_bases,
            source.use_basis,
            source.residual_scale_rationale,
        ),
        objective=ObjectiveContract(
            residual_family,
            "native_scaled_least_squares",
            "row scales retained by the canonical observations",
            "no covariance supplied",
            "squared",
            (),
            "fail",
        ),
        row_provenance={
            row.row_id: RowProvenance(
                AcquisitionClass.DIRECT_MEASUREMENT,
                "unique retained row",
                "included",
                "declared by source workflow",
                "not censored",
                "retained under the declared row policy",
            )
            for row in problem.observations
        },
    )
    return prepare_fit(
        model,
        datasets=(dataset,),
        parameters=tuple(
            ParameterRequest(
                coordinate.family,
                coordinate.identity,
                coordinate.transform,
                coordinate.lower_bound,
                coordinate.upper_bound,
            )
            for coordinate in problem.parameters
        ),
        parameter_slot_indices=problem.parameter_slot_indices,
        start_vectors=problem.start_vectors,
        solver=SolverControls(
            problem.maximum_iterations,
            problem.maximum_solver_time_seconds,
            problem.function_tolerance,
            problem.gradient_tolerance,
            problem.parameter_tolerance,
        ),
        rank=RankControls(problem.maximum_condition_number),
        confirmation=ConfirmationControls(
            problem.confirmation_parameter_scaled_max_delta,
            problem.confirmation_cost_relative_delta,
        ),
    )


def _row(
    row_id: str = "may2015-ch4-c2h6-002",
    partition: ObservationPartition = ObservationPartition.TRAINING,
) -> FixedCompositionVleObservation:
    return FixedCompositionVleObservation(
        row_id=row_id,
        source_id="may-2015",
        source_locator=f"evidence/may-2015-methane-ethane-vle.csv:{row_id}",
        component_ids=("methane", "ethane"),
        temperature_k=203.22,
        pressure_pa=2_124_000.0,
        liquid_mole_fraction_first=0.3653,
        vapor_mole_fraction_first=0.8667,
        pressure_scale_pa=2_124_000.0,
        chemical_potential_scales=(1.0, 1.0),
        liquid_volume_origin_m3_per_mol=6.0e-5,
        liquid_volume_start_m3_per_mol=6.5e-5,
        liquid_volume_bounds_m3_per_mol=(2.0e-5, 1.0e-4),
        vapor_volume_origin_m3_per_mol=9.0e-4,
        vapor_volume_start_m3_per_mol=1.0e-3,
        vapor_volume_bounds_m3_per_mol=(1.0e-4, 1.0e-2),
        partition=partition,
    )


def _problem(
    model: Mixture,
    rows: tuple[FixedCompositionVleObservation, ...] | None = None,
    family: ParameterFamily = ParameterFamily.K_IJ,
) -> RegressionProblem:
    capability = _capability(model, family)
    observations = rows or (_row(),)
    source = SourceDescriptor(
        source_id="may-2015",
        citation="May et al. (2015) methane/ethane VLE.",
        durable_locator="evidence/may-2015-methane-ethane-vle.csv",
        source_artifact_sha256=(
            "5cd1e74925a3c6504f5106dcf911f2cae2d6e99a5133fccc20454d8991bdbc7f"
        ),
        canonical_dataset_sha256=canonical_dataset_sha256(observations),
        transformation_record="No transformation for this contract test.",
        units_and_bases="T/K, P/Pa, mole fractions.",
        use_basis="Regression derivative contract evidence.",
        residual_scale_rationale="Pressure by observed P; mu/RT dimensionless.",
    )
    parameter = ParameterCoordinate(
        family=family,
        identity=PairParameterIdentity("methane", "ethane"),
        capability_id=capability.capability_id,
        provider_parameter_fingerprint=capability.parameter_fingerprint,
        provider_topology_fingerprint=capability.topology_fingerprint,
        unit="1",
        transform=AffineParameterTransform(origin=0.0, scale=0.01),
        lower_bound=-0.15,
        upper_bound=0.10,
    )
    return RegressionProblem(
        sources=(source,),
        parameters=(parameter,),
        parameter_slot_indices=(0,),
        start_vectors=((0.0,), (-0.05,), (0.05,)),
        observations=observations,
        maximum_condition_number=1.0e10,
        maximum_iterations=50,
        maximum_solver_time_seconds=30.0,
        function_tolerance=1.0e-12,
        gradient_tolerance=1.0e-12,
        parameter_tolerance=1.0e-12,
        confirmation_parameter_scaled_max_delta=1.0e-5,
        confirmation_cost_relative_delta=1.0e-8,
    )


def _mock_general_native_result() -> tuple[object, ...]:
    """Return one complete native result for status/diagnostic seam tests."""
    row = (
        "may2015-ch4-c2h6-002",
        "training",
        6.0e-5,
        9.0e-4,
        (0.0, 0.0, 0.0, 0.0),
        True,
        "",
    )
    return (
        "CONVERGENCE",
        True,
        1.0,
        0.1,
        2,
        -0.01,
        0.09,
        "",
        (0.0,) * 4,
        (0.0,) * 12,
        (2.0, 1.0, 0.5),
        3,
        2.0,
        (0.1,),
        1,
        1.0,
        2,
        0.0,
        0.0,
        True,
        (row,),
        "",
        3,
        4,
        4,
        3,
    )


def _general_result_signature(result: RegressionResult) -> tuple[object, ...]:
    return (
        result.parameters[0].final,
        result.final_cost,
        result.jacobian.residual_count,
        result.jacobian.variable_count,
        result.jacobian.full_rank,
        result.jacobian.projected_parameter_rank,
        result.training_row_count,
        result.held_out_row_count,
        result.stress_row_count,
        result.evaluated_row_count,
        result.skipped_row_count,
        result.failed_row_count,
    )


def _pure_problem(
    model: Mixture,
    family: ParameterFamily,
    *,
    all_training_rows: bool = False,
) -> RegressionProblem:
    dataset = load_pure_saturation_dataset("methane")
    source_rows = (
        dataset.training_rows if all_training_rows else dataset.training_rows[:1]
    )
    observations = tuple(
        PureSaturationObservation(
            row_id=row.row_id,
            source_id=row.source_id,
            source_locator=f"{dataset.source.locator}:{row.row_id}",
            component_id="methane",
            temperature_k=row.temperature_k,
            pressure_pa=row.pressure_pa,
            liquid_density_kg_per_m3=row.liquid_density_kg_m3,
            molar_mass_kg_per_mol=0.016043,
            pressure_scale_pa=2.0 * row.pressure_pa,
            chemical_potential_scale=2.0,
            liquid_density_scale_kg_per_m3=(2.0 * row.liquid_density_kg_m3),
            liquid_volume_origin_m3_per_mol=(0.016043 / row.liquid_density_kg_m3),
            liquid_volume_start_m3_per_mol=(0.016043 / row.liquid_density_kg_m3),
            liquid_volume_bounds_m3_per_mol=(2.0e-5, 1.0e-4),
            vapor_volume_origin_m3_per_mol=(
                8.31446261815324 * row.temperature_k / row.pressure_pa
            ),
            vapor_volume_start_m3_per_mol=(
                8.31446261815324 * row.temperature_k / row.pressure_pa
            ),
            vapor_volume_bounds_m3_per_mol=(1.5e-4, 0.1),
            partition=ObservationPartition.TRAINING,
        )
        for row in source_rows
    )
    capability = _capability(model, family)
    settings = {
        ParameterFamily.SEGMENT_COUNT: (1.0, 0.1, 0.5, 3.5, (1.0, 1.1), "1"),
        ParameterFamily.SEGMENT_DIAMETER: (
            3.7039,
            0.1,
            2.0,
            5.0,
            (3.7039, 3.6),
            "angstrom",
        ),
        ParameterFamily.DISPERSION_ENERGY_OVER_K: (
            150.03,
            10.0,
            50.0,
            400.0,
            (150.03, 160.0),
            "K",
        ),
    }
    origin, scale, lower, upper, starts, unit = settings[family]
    source = SourceDescriptor(
        source_id=dataset.source.source_id,
        citation=dataset.source.citation,
        durable_locator=dataset.source.locator,
        source_artifact_sha256=dataset.source.data_sha256,
        canonical_dataset_sha256=canonical_dataset_sha256(observations),
        transformation_record=dataset.source.transformation,
        units_and_bases="T/K, P/Pa, saturated liquid density/(kg/m3).",
        use_basis=dataset.source.use_basis,
        residual_scale_rationale=(
            "Preserve the accepted equal four-residual weighting: scales are "
            "2P, 2 for mu/RT, and twice observed liquid density."
        ),
    )
    parameter = ParameterCoordinate(
        family=family,
        identity=ComponentParameterIdentity("methane"),
        capability_id=capability.capability_id,
        provider_parameter_fingerprint=capability.parameter_fingerprint,
        provider_topology_fingerprint=capability.topology_fingerprint,
        unit=unit,
        transform=AffineParameterTransform(origin=origin, scale=scale),
        lower_bound=lower,
        upper_bound=upper,
    )
    return RegressionProblem(
        sources=(source,),
        parameters=(parameter,),
        parameter_slot_indices=(0,),
        start_vectors=tuple((start,) for start in starts),
        observations=observations,
        maximum_condition_number=1.0e12,
        maximum_iterations=500,
        maximum_solver_time_seconds=30.0,
        function_tolerance=1.0e-10,
        gradient_tolerance=1.0e-10,
        parameter_tolerance=1.0e-10,
        confirmation_parameter_scaled_max_delta=1.0e-5,
        confirmation_cost_relative_delta=1.0e-8,
    )


def _joint_pure_problem(model: Mixture) -> RegressionProblem:
    families = (
        ParameterFamily.SEGMENT_COUNT,
        ParameterFamily.SEGMENT_DIAMETER,
        ParameterFamily.DISPERSION_ENERGY_OVER_K,
    )
    scalar_problems = tuple(
        _pure_problem(model, family, all_training_rows=True) for family in families
    )
    return replace(
        scalar_problems[0],
        parameters=tuple(problem.parameters[0] for problem in scalar_problems),
        parameter_slot_indices=(0, 1, 2),
        start_vectors=(
            (1.08, 3.555744, 157.5315),
            (1.0, 3.7, 150.0),
        ),
    )


def _solvation_factor_problem(model: Mixture) -> RegressionProblem:
    specification = FIGIEL_WATER_SOLVATION_FACTOR_V1
    capability = _capability(model, ParameterFamily.SOLVATION_FACTOR)
    observations = tuple(
        MeanIonicActivityObservation(
            row_id=row.row_id,
            source_id="hamer-wu-1972-nabr",
            source_locator=(
                f"validation:data/figiel-2025-hamer-wu-miac.csv:{row.row_id}"
            ),
            component_ids=capability.component_ids,
            active_component_id="water",
            temperature_k=specification.temperature_k,
            pressure_pa=specification.pressure_pa,
            formula_unit_molality_mol_per_kg=(row.molality_mol_per_kg),
            observed_mean_ionic_activity_coefficient=row.gamma_pm_m,
            relative_residual_scale=1.0,
            partition=ObservationPartition.TRAINING,
        )
        for row in specification.observations
    )
    source = SourceDescriptor(
        source_id="hamer-wu-1972-nabr",
        citation=(
            "Hamer and Wu (1972), NaBr mean ionic activity coefficients; "
            "audited by Validation packet 8944d34f."
        ),
        durable_locator=("validation:data/figiel-2025-hamer-wu-miac.csv"),
        source_artifact_sha256=specification.source_hamer_wu_csv_sha256,
        canonical_dataset_sha256=canonical_dataset_sha256(observations),
        transformation_record=(
            "Selected all 21 audited NaBr rows; no row duplication or "
            "response transformation."
        ),
        units_and_bases=(
            "T/K, P/Pa, formula-unit molality/(mol/kg), gamma_pm on the molality basis."
        ),
        use_basis="All rows are in-sample parameter-recovery targets.",
        residual_scale_rationale=(
            "Preserve the frozen dimensionless residual 1 - gamma_model/gamma_observed."
        ),
    )
    parameter = ParameterCoordinate(
        family=ParameterFamily.SOLVATION_FACTOR,
        identity=ComponentParameterIdentity("water"),
        capability_id=capability.capability_id,
        provider_parameter_fingerprint=capability.parameter_fingerprint,
        provider_topology_fingerprint=capability.topology_fingerprint,
        unit="1",
        transform=AffineParameterTransform(origin=1.5, scale=0.1),
        lower_bound=specification.parameter_bounds[0],
        upper_bound=specification.parameter_bounds[1],
    )
    return RegressionProblem(
        sources=(source,),
        parameters=(parameter,),
        parameter_slot_indices=(0,),
        start_vectors=tuple((start,) for start in specification.starts),
        observations=observations,
        maximum_condition_number=1.0e10,
        maximum_iterations=specification.max_num_iterations,
        maximum_solver_time_seconds=180.0,
        function_tolerance=specification.function_tolerance,
        gradient_tolerance=specification.gradient_tolerance,
        parameter_tolerance=specification.parameter_tolerance,
        confirmation_parameter_scaled_max_delta=1.0e-4,
        confirmation_cost_relative_delta=1.0e-8,
    )


def _aqueous_kij_problem(model: Mixture) -> RegressionProblem:
    specification = FIGIEL_AQUEOUS_KIJ_V1
    capability = next(
        capability
        for capability in parameter_capabilities(model)
        if not isinstance(capability, UnsupportedParameterCapability)
        and capability.capability_id == "aqueous_water_cation_kij_miac_v1"
    )
    observations = tuple(
        AqueousKijMeanIonicActivityObservation(
            row_id=row.row_id,
            source_id="hamer-wu-1972-nabr",
            source_locator=(
                f"validation:data/hamer-wu-1972-aqueous-alkali-halides.csv:{row.row_id}"
            ),
            component_ids=capability.component_ids,
            active_pair_component_ids=("water", "sodium-cation"),
            fixed_k_ij=(-0.3, -0.3, 0.65),
            temperature_k=specification.temperature_k,
            pressure_pa=specification.pressure_pa,
            formula_unit_molality_mol_per_kg=row.molality_mol_per_kg,
            observed_mean_ionic_activity_coefficient=row.gamma_pm_m,
            relative_residual_scale=1.0,
            partition=ObservationPartition.TRAINING,
        )
        for row in specification.observations
        if row.salt == "NaBr"
    )
    source = SourceDescriptor(
        source_id="hamer-wu-1972-nabr",
        citation=(
            "Hamer and Wu (1972), NaBr mean ionic activity coefficients; "
            "audited by Validation packet 8944d34f."
        ),
        durable_locator=("validation:data/hamer-wu-1972-aqueous-alkali-halides.csv"),
        source_artifact_sha256=specification.source_hamer_wu_csv_sha256,
        canonical_dataset_sha256=canonical_dataset_sha256(observations),
        transformation_record=(
            "Selected all 21 audited NaBr rows; fixed water-anion and "
            "cation-anion interactions are explicit workflow inputs."
        ),
        units_and_bases=(
            "T/K, P/Pa, formula-unit molality/(mol/kg), gamma_pm on the "
            "molality basis; k_ij dimensionless."
        ),
        use_basis="All rows are in-sample scalar parameter-recovery targets.",
        residual_scale_rationale=(
            "Use the frozen dimensionless residual 1 - gamma_model/gamma_observed."
        ),
    )
    parameter = ParameterCoordinate(
        family=ParameterFamily.K_IJ,
        identity=PairParameterIdentity("water", "sodium-cation"),
        capability_id=capability.capability_id,
        provider_parameter_fingerprint=capability.parameter_fingerprint,
        provider_topology_fingerprint=capability.topology_fingerprint,
        unit="1",
        transform=AffineParameterTransform(origin=-0.3, scale=0.1),
        lower_bound=-1.0,
        upper_bound=1.0,
    )
    return RegressionProblem(
        sources=(source,),
        parameters=(parameter,),
        parameter_slot_indices=(0,),
        start_vectors=((0.0,), (0.25,)),
        observations=observations,
        maximum_condition_number=1.0e10,
        maximum_iterations=50,
        maximum_solver_time_seconds=180.0,
        function_tolerance=1.0e-10,
        gradient_tolerance=1.0e-10,
        parameter_tolerance=1.0e-10,
        confirmation_parameter_scaled_max_delta=1.0e-5,
        confirmation_cost_relative_delta=1.0e-8,
    )


def _born_diameter_problem(model: Mixture, target_index: int) -> RegressionProblem:
    specification = FIGIEL_BORN_DIAMETER_TRACER_V1
    target = specification.targets[target_index]
    capability = _capability(model, ParameterFamily.BORN_DIAMETER)
    observation = SolvationGibbsObservation(
        row_id=target.target_id,
        source_id="figiel-2025-s5-reported-averages",
        source_locator=(f"{specification.source_locator}:{target.target_id}"),
        component_ids=target.component_order,
        active_component_id=target.active_component_id,
        temperature_k=specification.temperature_k,
        pressure_pa=specification.pressure_pa,
        observed_solvation_gibbs_j_per_mol=target.target_j_per_mol,
        residual_scale_j_per_mol=abs(target.target_j_per_mol),
        partition=ObservationPartition.TRAINING,
    )
    source = SourceDescriptor(
        source_id="figiel-2025-s5-reported-averages",
        citation=(
            "Figiel, Yu, and Held (2025), SI Table S5 reported-average "
            "water solvation Gibbs energies."
        ),
        durable_locator=specification.source_si_doi,
        source_artifact_sha256=specification.packaged_targets_sha256,
        canonical_dataset_sha256=canonical_dataset_sha256((observation,)),
        transformation_record=(
            "Converted the selected reported average from kJ/mol to J/mol; "
            "the 27 supporting rows are provenance, not residual duplication."
        ),
        units_and_bases=(
            f"T/K, P/Pa, solvation Gibbs energy/(J/mol); {specification.source_basis}."
        ),
        use_basis=("One in-sample source-bound Born-diameter recovery problem."),
        residual_scale_rationale=(
            "Scale the Gibbs-energy difference by the magnitude of its "
            "reported target, preserving the frozen tracer contract."
        ),
    )
    parameter = ParameterCoordinate(
        family=ParameterFamily.BORN_DIAMETER,
        identity=ComponentParameterIdentity(target.active_component_id),
        capability_id=capability.capability_id,
        provider_parameter_fingerprint=capability.parameter_fingerprint,
        provider_topology_fingerprint=capability.topology_fingerprint,
        unit="angstrom",
        transform=AffineParameterTransform(
            origin=specification.diameter_origin_angstrom,
            scale=specification.diameter_scale_angstrom,
        ),
        lower_bound=specification.diameter_bounds_angstrom[0],
        upper_bound=specification.diameter_bounds_angstrom[1],
    )
    return RegressionProblem(
        sources=(source,),
        parameters=(parameter,),
        parameter_slot_indices=(0,),
        start_vectors=tuple(
            (start[target_index],) for start in specification.start_diameters_angstrom
        ),
        observations=(observation,),
        maximum_condition_number=1.0e10,
        maximum_iterations=specification.max_num_iterations,
        maximum_solver_time_seconds=180.0,
        function_tolerance=specification.function_tolerance,
        gradient_tolerance=specification.gradient_tolerance,
        parameter_tolerance=specification.parameter_tolerance,
        confirmation_parameter_scaled_max_delta=(
            specification.confirmation_parameter_scaled_max_delta
        ),
        confirmation_cost_relative_delta=1.0e-8,
    )


def _ionic_region_permittivity_problem(model: Mixture) -> RegressionProblem:
    specification = FIGIEL_BORN_DIAMETER_TRACER_V1
    target = specification.targets[1]
    capability = _capability(model, ParameterFamily.IONIC_REGION_RELATIVE_PERMITTIVITY)
    observation = SolvationGibbsObservation(
        row_id=target.target_id,
        source_id="figiel-2025-s5-reported-averages",
        source_locator=f"{specification.source_locator}:{target.target_id}",
        component_ids=target.component_order,
        active_component_id=target.active_component_id,
        temperature_k=specification.temperature_k,
        pressure_pa=specification.pressure_pa,
        observed_solvation_gibbs_j_per_mol=target.target_j_per_mol,
        residual_scale_j_per_mol=abs(target.target_j_per_mol),
        partition=ObservationPartition.TRAINING,
    )
    source = SourceDescriptor(
        source_id="figiel-2025-s5-reported-averages",
        citation=(
            "Figiel, Yu, and Held (2025), SI Table S5 reported-average "
            "water solvation Gibbs energies."
        ),
        durable_locator=specification.source_si_doi,
        source_artifact_sha256=specification.packaged_targets_sha256,
        canonical_dataset_sha256=canonical_dataset_sha256((observation,)),
        transformation_record=(
            "Converted the sodium reported average from kJ/mol to J/mol; "
            "the 27 supporting rows are provenance, not residual duplication."
        ),
        units_and_bases=(
            "T/K, P/Pa, solvation Gibbs energy/(J/mol), and dimensionless "
            "ionic-region relative permittivity."
        ),
        use_basis=(
            "One in-sample derivative and scalar-fit capability check; "
            "Figiel reports the permittivity as fixed model input."
        ),
        residual_scale_rationale=(
            "Scale the Gibbs-energy difference by the reported-target "
            "magnitude; this is a numerical scale, not uncertainty."
        ),
    )
    parameter = ParameterCoordinate(
        family=ParameterFamily.IONIC_REGION_RELATIVE_PERMITTIVITY,
        identity=ModelParameterIdentity(),
        capability_id=capability.capability_id,
        provider_parameter_fingerprint=capability.parameter_fingerprint,
        provider_topology_fingerprint=capability.topology_fingerprint,
        unit="1",
        transform=AffineParameterTransform(origin=8.0, scale=2.0),
        lower_bound=1.01,
        upper_bound=50.0,
    )
    return RegressionProblem(
        sources=(source,),
        parameters=(parameter,),
        parameter_slot_indices=(0,),
        start_vectors=((4.0,), (12.0,)),
        observations=(observation,),
        maximum_condition_number=1.0e10,
        maximum_iterations=50,
        maximum_solver_time_seconds=180.0,
        function_tolerance=1.0e-12,
        gradient_tolerance=1.0e-12,
        parameter_tolerance=1.0e-12,
        confirmation_parameter_scaled_max_delta=1.0e-5,
        confirmation_cost_relative_delta=1.0e-8,
    )


def _solvent_relative_permittivity_problem(model: Mixture) -> RegressionProblem:
    specification = FIGIEL_BORN_DIAMETER_TRACER_V1
    target = specification.targets[1]
    capability = _capability(model, ParameterFamily.RELATIVE_PERMITTIVITY)
    observation = SolvationGibbsObservation(
        row_id=target.target_id,
        source_id="figiel-2025-s5-reported-averages",
        source_locator=f"{specification.source_locator}:{target.target_id}",
        component_ids=target.component_order,
        active_component_id=target.active_component_id,
        temperature_k=specification.temperature_k,
        pressure_pa=specification.pressure_pa,
        observed_solvation_gibbs_j_per_mol=target.target_j_per_mol,
        residual_scale_j_per_mol=abs(target.target_j_per_mol),
        partition=ObservationPartition.TRAINING,
    )
    source = SourceDescriptor(
        source_id="figiel-2025-s5-reported-averages",
        citation=(
            "Figiel, Yu, and Held (2025), SI Table S5 reported-average "
            "water solvation Gibbs energies."
        ),
        durable_locator=specification.source_si_doi,
        source_artifact_sha256=specification.packaged_targets_sha256,
        canonical_dataset_sha256=canonical_dataset_sha256((observation,)),
        transformation_record=(
            "Converted the sodium reported average from kJ/mol to J/mol; "
            "the 27 supporting rows are provenance, not residual duplication."
        ),
        units_and_bases=(
            "T/K, P/Pa, solvation Gibbs energy/(J/mol), and dimensionless "
            "water relative permittivity."
        ),
        use_basis=(
            "One in-sample derivative and scalar-fit capability check; "
            "Figiel reports the water permittivity as fixed model input."
        ),
        residual_scale_rationale=(
            "Scale the Gibbs-energy difference by the reported-target "
            "magnitude; this is a numerical scale, not uncertainty."
        ),
    )
    parameter = ParameterCoordinate(
        family=ParameterFamily.RELATIVE_PERMITTIVITY,
        identity=ComponentParameterIdentity("water"),
        capability_id=capability.capability_id,
        provider_parameter_fingerprint=capability.parameter_fingerprint,
        provider_topology_fingerprint=capability.topology_fingerprint,
        unit="1",
        transform=AffineParameterTransform(origin=78.09, scale=10.0),
        lower_bound=1.01,
        upper_bound=200.0,
    )
    return RegressionProblem(
        sources=(source,),
        parameters=(parameter,),
        parameter_slot_indices=(0,),
        start_vectors=((50.0,), (110.0,)),
        observations=(observation,),
        maximum_condition_number=1.0e10,
        maximum_iterations=50,
        maximum_solver_time_seconds=180.0,
        function_tolerance=1.0e-12,
        gradient_tolerance=1.0e-12,
        parameter_tolerance=1.0e-12,
        confirmation_parameter_scaled_max_delta=1.0e-5,
        confirmation_cost_relative_delta=1.0e-8,
    )


def _dielectric_suppression_problem(model: Mixture) -> RegressionProblem:
    capability = _capability(
        model,
        ParameterFamily.ION_FRACTION_SUPPRESSION_COEFFICIENT,
    )
    source_rows = (
        ("figiel2025-dielectric-water-008", 0.010880829, 71.61417323),
        ("figiel2025-dielectric-water-009", 0.022927461, 66.65354331),
        ("figiel2025-dielectric-water-010", 0.033808290, 62.87401575),
    )
    observations = tuple(
        RelativePermittivityRatioObservation(
            row_id=row_id,
            source_id="retained-lab-figiel-figure-2",
            source_locator=f"validation:figiel-ledger:{row_id}",
            solvent_id="water",
            component_ids=capability.component_ids,
            temperature_k=298.15,
            pressure_pa=100_000.0,
            total_ion_mole_fraction=ion_fraction,
            observed_relative_permittivity_ratio=permittivity / 78.09,
            residual_scale=1.0,
            partition=ObservationPartition.TRAINING,
        )
        for row_id, ion_fraction, permittivity in source_rows
    )
    source = SourceDescriptor(
        source_id="retained-lab-figiel-figure-2",
        citation=(
            "Figiel, Yu, and Held (2025), Figure 2 water dielectric "
            "correlation; retained digitization."
        ),
        durable_locator="validation:data/figiel-2025-regression-target-ledger.csv",
        source_artifact_sha256=(
            "09e1e820a55861b835fdab27df5134451ecc9329c6d512dcf26565b267b387a6"
        ),
        canonical_dataset_sha256=canonical_dataset_sha256(observations),
        transformation_record=(
            "Selected three NaBr rows without response transformation."
        ),
        units_and_bases=(
            "T/K, P/Pa, total ion mole fraction, dimensionless static "
            "relative dielectric constant."
        ),
        use_basis="Digitized in-sample implementation evidence only.",
        residual_scale_rationale=(
            "Equal raw dimensionless relative-permittivity residuals; no "
            "pointwise uncertainty is available."
        ),
    )
    parameter = ParameterCoordinate(
        family=ParameterFamily.ION_FRACTION_SUPPRESSION_COEFFICIENT,
        identity=ModelParameterIdentity(),
        capability_id=capability.capability_id,
        provider_parameter_fingerprint=capability.parameter_fingerprint,
        provider_topology_fingerprint=capability.topology_fingerprint,
        unit="1",
        transform=AffineParameterTransform(origin=7.0, scale=1.0),
        lower_bound=0.01,
        upper_bound=30.0,
    )
    return RegressionProblem(
        sources=(source,),
        parameters=(parameter,),
        parameter_slot_indices=(0,),
        start_vectors=((2.0,), (12.0,)),
        observations=observations,
        maximum_condition_number=1.0e10,
        maximum_iterations=50,
        maximum_solver_time_seconds=180.0,
        function_tolerance=1.0e-12,
        gradient_tolerance=1.0e-12,
        parameter_tolerance=1.0e-12,
        confirmation_parameter_scaled_max_delta=1.0e-5,
        confirmation_cost_relative_delta=1.0e-8,
    )


def _ion_solvation_kij_problem(
    model: Mixture,
    *,
    capability_id: str = "ion_solvation_solvent_cation_kij_v1",
    active_component_id: str = "potassium-cation",
    active_pair_component_ids: tuple[str, str] = (
        "methanol",
        "potassium-cation",
    ),
) -> RegressionProblem:
    capability = next(
        capability
        for capability in parameter_capabilities(model)
        if not isinstance(capability, UnsupportedParameterCapability)
        and capability.capability_id == capability_id
    )
    targets = (("figiel2025-constructed-gsolv-Kp-methanol-011", -298.25858),)
    observations = tuple(
        IonSolvationKijObservation(
            row_id=row_id,
            source_id="figiel-constructed-k-methanol",
            source_locator=f"validation:figiel-ledger:{row_id}",
            component_ids=capability.component_ids,
            active_component_id=active_component_id,
            active_pair_component_ids=active_pair_component_ids,
            fixed_k_ij=(0.32, 0.15, -0.35),
            temperature_k=298.15,
            pressure_pa=100_000.0,
            observed_solvation_gibbs_j_per_mol=value * 1000.0,
            residual_scale_j_per_mol=300_000.0,
            partition=ObservationPartition.TRAINING,
        )
        for row_id, value in targets
    )
    source = SourceDescriptor(
        source_id="figiel-constructed-k-methanol",
        citation="Figiel, Yu, and Held (2025), equation 19 and Figure 6.",
        durable_locator="validation:data/figiel-2025-regression-target-ledger.csv",
        source_artifact_sha256=(
            "f405a3e48d21cd979a8dd480d5f8cb3be40754f5d6babf368b505b5f305607f0"
        ),
        canonical_dataset_sha256=canonical_dataset_sha256(observations),
        transformation_record=(
            "Converted the retained constructed kJ/mol targets to J/mol."
        ),
        units_and_bases="T/K, P/Pa, x-process ion solvation Gibbs/J/mol.",
        use_basis=(
            "Nearest digitized pure-methanol endpoint used as constructed "
            "in-sample implementation evidence only."
        ),
        residual_scale_rationale=(
            "A declared common 300 kJ/mol magnitude scale; no pointwise "
            "uncertainty is available."
        ),
    )
    parameter = ParameterCoordinate(
        family=ParameterFamily.K_IJ,
        identity=PairParameterIdentity(*active_pair_component_ids),
        capability_id=capability.capability_id,
        provider_parameter_fingerprint=capability.parameter_fingerprint,
        provider_topology_fingerprint=capability.topology_fingerprint,
        unit="1",
        transform=AffineParameterTransform(origin=0.3, scale=0.1),
        lower_bound=-1.5,
        upper_bound=1.5,
    )
    return RegressionProblem(
        sources=(source,),
        parameters=(parameter,),
        parameter_slot_indices=(0,),
        start_vectors=((0.0,), (0.7,)),
        observations=observations,
        maximum_condition_number=1.0e10,
        maximum_iterations=50,
        maximum_solver_time_seconds=180.0,
        function_tolerance=1.0e-12,
        gradient_tolerance=1.0e-12,
        parameter_tolerance=1.0e-12,
        confirmation_parameter_scaled_max_delta=1.0e-5,
        confirmation_cost_relative_delta=1.0e-8,
    )


def _replace_observations(
    problem: RegressionProblem, observations: tuple[object, ...]
) -> RegressionProblem:
    source = replace(
        problem.sources[0],
        canonical_dataset_sha256=canonical_dataset_sha256(observations),
    )
    return replace(problem, sources=(source,), observations=observations)


def test_installed_provider_advertises_exact_neutral_binary_kij_contract(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    model = _model()

    capability = _capability(model)

    assert capability.capability_id == "neutral_binary_phase_kij_v1"
    assert capability.family is ParameterFamily.K_IJ
    assert capability.component_ids == ("methane", "ethane")
    assert capability.coordinate_kinds == ("amount", "amount", "volume", "k_ij")
    assert capability.coordinate_units == ("mol", "mol", "m3", "dimensionless")
    assert capability.parameter_fingerprint == model.parameter_fingerprint
    assert capability.topology_fingerprint.startswith("sha256:")
    assert capability.temperature_min_k <= 203.22 <= capability.temperature_max_k
    assert capability.derivative_order == 2
    assert capability.maturity == "DERIVATIVE_READY"
    assert capability.authority_effect == "NONE"
    assert capability.identity_shape == "unordered_component_pair"
    assert capability.observation_contract == "fixed_composition_helmholtz_phase"
    assert capability.model_domain == "neutral_nonassociating_binary"
    assert capability.tensor_layout == "row_major"
    assert capability.state_coordinate_count == 3
    assert capability.active_parameter_count == 1
    assert capability.helmholtz_basis_id == (
        "A_over_RT_reference_amount:n_ref=1mol:rho_ref=1mol_per_m3"
    )
    assert capability.unsupported_status == "UNSUPPORTED_MODEL"
    assert capability.domain_status == "DOMAIN_ERROR"

    lij = _capability(model, ParameterFamily.L_IJ)
    assert lij.capability_id == "neutral_binary_phase_lij_v1"
    assert lij.component_ids == capability.component_ids
    assert lij.coordinate_kinds == ("amount", "amount", "volume", "l_ij")
    assert lij.coordinate_units == capability.coordinate_units
    assert lij.parameter_fingerprint == capability.parameter_fingerprint
    assert lij.topology_fingerprint == capability.topology_fingerprint
    assert lij.derivative_order == 2
    assert lij.observation_contract == capability.observation_contract
    assert lij.model_domain == capability.model_domain

    raw = parameter_regression._native.parameter_capabilities(native_sdk(model))
    monkeypatch.setattr(
        parameter_regression._native,
        "parameter_capabilities",
        lambda _: (*raw, (999, 1, 998)),
    )
    known = parameter_capabilities(model)
    assert known[0] == capability
    unsupported = known[-1]
    assert unsupported == UnsupportedParameterCapability(
        capability_code=999,
        schema_version=1,
        parameter_family_code=998,
    )


def test_installed_provider_advertises_exact_direct_observable_contracts() -> None:
    solvation = _capability(_aqueous_model(), ParameterFamily.SOLVATION_FACTOR)
    born = _capability(
        _aqueous_model(("water", "sodium-cation", "chloride-anion")),
        ParameterFamily.BORN_DIAMETER,
    )

    assert (
        solvation.capability_id,
        solvation.component_ids,
        solvation.coordinate_kinds,
        solvation.coordinate_units,
        solvation.derivative_order,
        solvation.identity_shape,
        solvation.observation_contract,
        solvation.model_domain,
        solvation.state_coordinate_count,
        solvation.active_parameter_count,
        solvation.active_component_ids,
    ) == (
        "aqueous_solvation_factor_miac_v1",
        ("water", "sodium-cation", "bromide-anion"),
        ("solvation_factor",),
        ("dimensionless",),
        1,
        "component",
        "aqueous_mean_ionic_activity",
        "figiel_aqueous_nabr",
        0,
        1,
        ("water",),
    )
    assert (
        born.capability_id,
        born.coordinate_kinds,
        born.coordinate_units,
        born.derivative_order,
        born.observation_contract,
        born.model_domain,
        born.active_component_ids,
    ) == (
        "ion_solvation_born_v1",
        ("born_diameter",),
        ("angstrom",),
        1,
        "ion_solvation_gibbs",
        "figiel_water_single_ion",
        ("sodium-cation",),
    )

    dielectric = _capability(
        _aqueous_model(),
        ParameterFamily.ION_FRACTION_SUPPRESSION_COEFFICIENT,
    )
    assert (
        dielectric.capability_id,
        dielectric.coordinate_kinds,
        dielectric.coordinate_units,
        dielectric.derivative_order,
        dielectric.identity_shape,
        dielectric.observation_contract,
        dielectric.model_domain,
        dielectric.active_component_ids,
    ) == (
        "ion_fraction_suppression_v1",
        ("ion_fraction_suppression_coefficient",),
        ("dimensionless",),
        1,
        "model",
        "relative_permittivity_ratio",
        "figiel_dielectric",
        (),
    )

    ionic_permittivity = _capability(
        _aqueous_model(("water", "sodium-cation", "chloride-anion")),
        ParameterFamily.IONIC_REGION_RELATIVE_PERMITTIVITY,
    )
    assert (
        ionic_permittivity.capability_id,
        ionic_permittivity.coordinate_kinds,
        ionic_permittivity.coordinate_units,
        ionic_permittivity.derivative_order,
        ionic_permittivity.identity_shape,
        ionic_permittivity.observation_contract,
        ionic_permittivity.model_domain,
        ionic_permittivity.active_component_ids,
    ) == (
        "ion_solvation_ionic_region_permittivity_v1",
        ("ionic_region_relative_permittivity",),
        ("dimensionless",),
        1,
        "model",
        "ion_solvation_gibbs",
        "figiel_water_single_ion",
        (),
    )

    solvent_permittivity = _capability(
        _aqueous_model(("water", "sodium-cation", "chloride-anion")),
        ParameterFamily.RELATIVE_PERMITTIVITY,
    )
    assert (
        solvent_permittivity.capability_id,
        solvent_permittivity.coordinate_kinds,
        solvent_permittivity.coordinate_units,
        solvent_permittivity.derivative_order,
        solvent_permittivity.identity_shape,
        solvent_permittivity.observation_contract,
        solvent_permittivity.model_domain,
        solvent_permittivity.active_component_ids,
    ) == (
        "ion_solvation_solvent_permittivity_v1",
        ("relative_permittivity",),
        ("dimensionless",),
        1,
        "component",
        "ion_solvation_gibbs",
        "figiel_water_single_ion",
        ("water",),
    )


def test_installed_provider_advertises_each_aqueous_kij_miac_contract() -> None:
    model = _aqueous_kij_models(FIGIEL_AQUEOUS_KIJ_V1)[4]
    capabilities = tuple(
        capability
        for capability in parameter_capabilities(model)
        if not isinstance(capability, UnsupportedParameterCapability)
        and capability.family is ParameterFamily.K_IJ
        and capability.observation_contract == "aqueous_mean_ionic_activity"
    )

    assert tuple(capability.capability_id for capability in capabilities) == (
        "aqueous_water_cation_kij_miac_v1",
        "aqueous_water_anion_kij_miac_v1",
        "aqueous_cation_anion_kij_miac_v1",
    )
    assert tuple(capability.active_component_ids for capability in capabilities) == (
        ("sodium-cation", "water"),
        ("bromide-anion", "water"),
        ("bromide-anion", "sodium-cation"),
    )
    assert all(
        capability.component_ids == ("water", "sodium-cation", "bromide-anion")
        and capability.coordinate_kinds == ("k_ij",)
        and capability.coordinate_units == ("dimensionless",)
        and capability.derivative_order == 1
        and capability.identity_shape == "unordered_component_pair"
        and capability.model_domain == "figiel_aqueous_nabr"
        for capability in capabilities
    )


@pytest.mark.campaign
def test_general_engine_fits_water_solvation_factor_over_all_nabr_rows() -> None:
    model = _fixed_water_factor_model(FIGIEL_WATER_SOLVATION_FACTOR_V1)
    result = fit_parameters(_solvation_factor_problem(model), model)

    assert result.solver_converged
    assert result.numerically_converged
    assert result.workflow_valid
    assert result.parameters[0].final == pytest.approx(
        1.5590515389548207, rel=1.0e-11, abs=1.0e-11
    )
    assert result.parameters[0].active_bound is None
    assert result.jacobian.residual_count == 21
    assert result.jacobian.variable_count == 1
    assert result.jacobian.full_rank == 1
    assert result.jacobian.projected_parameter_rank == 1
    assert result.confirmations_usable
    assert all(
        isinstance(row, DirectObservationRowDiagnostic)
        and row.evaluated
        and row.derivative_status == "EXACT_PROVIDER_FIRST_DERIVATIVE"
        for row in result.rows
    )


@pytest.mark.campaign
def test_general_engine_fits_organic_ion_solvation_kij_endpoint() -> None:
    model = _aqueous_model(("methanol", "potassium-cation", "bromide-anion"))
    result = fit_parameters(_ion_solvation_kij_problem(model), model)

    assert result.solver_converged
    assert result.numerically_converged
    assert result.workflow_valid
    assert result.parameters[0].final == pytest.approx(
        0.3467279724950645, rel=0.0, abs=2.0e-12
    )
    assert result.parameters[0].active_bound is None
    assert result.jacobian.full_rank == 1
    assert result.jacobian.projected_parameter_rank == 1
    assert result.confirmations_usable
    assert result.rows[0].derivative_status == ("EXACT_PROVIDER_FIRST_DERIVATIVE")


def test_general_engine_fits_dielectric_suppression_from_user_rows() -> None:
    model = _aqueous_model()
    result = fit_parameters(_dielectric_suppression_problem(model), model)

    assert result.solver_converged
    assert result.numerically_converged
    assert result.workflow_valid
    assert result.parameters[0].active_bound is None
    assert result.jacobian.residual_count == 3
    assert result.jacobian.variable_count == 1
    assert result.jacobian.full_rank == 1
    assert result.jacobian.projected_parameter_rank == 1
    assert result.confirmations_usable
    assert all(
        isinstance(row, DirectObservationRowDiagnostic)
        and row.observable == "relative_permittivity_ratio"
        and row.observable_unit == "1"
        and row.evaluated
        and row.derivative_status == "EXACT_PROVIDER_FIRST_DERIVATIVE"
        for row in result.rows
    )


@pytest.mark.campaign
def test_general_engine_fits_one_aqueous_kij_from_user_rows() -> None:
    model = _aqueous_kij_models(FIGIEL_AQUEOUS_KIJ_V1)[4]
    problem = _aqueous_kij_problem(model)
    observations_before = problem.observations
    fixed_context_before = tuple(row.fixed_k_ij for row in problem.observations)
    assert len(fixed_context_before) == 21

    result = fit_parameters(problem, model)

    assert result.solver_converged
    assert result.numerically_converged
    assert result.workflow_valid
    assert result.parameters[0].active_bound is None
    assert len(result.problem.parameters) == 1
    assert result.jacobian.residual_count == 21
    assert result.jacobian.variable_count == 1
    assert result.jacobian.full_rank == 1
    assert result.jacobian.projected_parameter_rank == 1
    assert result.confirmations_usable
    assert result.training_row_count == 21
    assert result.held_out_row_count == 0
    assert result.stress_row_count == 0
    assert result.evaluated_row_count == 21
    assert result.skipped_row_count == 0
    assert result.failed_row_count == 0
    assert problem.observations == observations_before
    assert tuple(row.fixed_k_ij for row in problem.observations) == fixed_context_before
    assert all(
        isinstance(row, DirectObservationRowDiagnostic)
        and row.evaluated
        and row.derivative_status == "EXACT_PROVIDER_FIRST_DERIVATIVE"
        for row in result.rows
    )


@pytest.mark.campaign
def test_ionic_region_permittivity_has_exact_rank_one_fit() -> None:
    target = FIGIEL_BORN_DIAMETER_TRACER_V1.targets[1]
    model = _aqueous_model(target.component_order)
    problem = _ionic_region_permittivity_problem(model)
    result = fit_parameters(problem, model)
    assert result.solver_converged
    assert result.numerically_converged
    assert result.workflow_valid
    assert result.parameters[0].active_bound is None
    assert 1.01 < result.parameters[0].final < 50.0
    assert result.jacobian.full_rank == 1
    assert result.jacobian.projected_parameter_rank == 1
    assert result.confirmations_usable
    assert isinstance(result.rows[0], DirectObservationRowDiagnostic)
    assert abs(result.rows[0].scaled_residual) <= 1.0e-10
    assert result.rows[0].derivative_status == ("EXACT_PROVIDER_FIRST_DERIVATIVE")


def test_ionic_region_permittivity_rejects_a_mislabeled_observable_ion() -> None:
    model = _aqueous_model(("water", "sodium-cation", "chloride-anion"))
    problem = _ionic_region_permittivity_problem(model)
    mislabeled = _replace_observations(
        problem,
        (
            replace(
                problem.observations[0],
                active_component_id="chloride-anion",
            ),
        ),
    )

    with pytest.raises(ValueError, match="active ion"):
        fit_parameters(mislabeled, model)


@pytest.mark.campaign
def test_solvent_relative_permittivity_has_exact_rank_one_fit() -> None:
    target = FIGIEL_BORN_DIAMETER_TRACER_V1.targets[1]
    model = _aqueous_model(target.component_order)
    problem = _solvent_relative_permittivity_problem(model)
    result = fit_parameters(problem, model)
    assert result.solver_converged
    assert result.numerically_converged
    assert result.workflow_valid
    assert result.parameters[0].active_bound is None
    assert 1.01 < result.parameters[0].final < 200.0
    assert result.jacobian.full_rank == 1
    assert result.jacobian.projected_parameter_rank == 1
    assert result.confirmations_usable
    assert isinstance(result.rows[0], DirectObservationRowDiagnostic)
    assert abs(result.rows[0].scaled_residual) <= 1.0e-10
    assert result.rows[0].derivative_status == ("EXACT_PROVIDER_FIRST_DERIVATIVE")


def test_solvent_relative_permittivity_rejects_a_mislabeled_observable_ion() -> None:
    model = _aqueous_model(("water", "sodium-cation", "chloride-anion"))
    problem = _solvent_relative_permittivity_problem(model)
    mislabeled = _replace_observations(
        problem,
        (
            replace(
                problem.observations[0],
                active_component_id="chloride-anion",
            ),
        ),
    )

    with pytest.raises(ValueError, match="active ion"):
        fit_parameters(mislabeled, model)


@pytest.mark.campaign
@pytest.mark.parametrize(
    ("target_index", "expected"),
    tuple(
        enumerate(
            (
                2.7888130173797934,
                3.4524616464076425,
                4.147266741279482,
                4.101505615791675,
                4.476998527506598,
            )
        )
    ),
)
def test_general_engine_fits_each_born_diameter_independently(
    target_index: int, expected: float
) -> None:
    target = FIGIEL_BORN_DIAMETER_TRACER_V1.targets[target_index]
    model = _aqueous_model(target.component_order)
    result = fit_parameters(_born_diameter_problem(model, target_index), model)

    assert result.solver_converged
    assert result.numerically_converged
    assert result.workflow_valid
    assert result.parameters[0].final == pytest.approx(
        expected, rel=2.0e-11, abs=2.0e-11
    )
    assert result.parameters[0].active_bound is None
    assert result.jacobian.residual_count == 1
    assert result.jacobian.variable_count == 1
    assert result.jacobian.full_rank == 1
    assert result.jacobian.projected_parameter_rank == 1
    assert result.confirmations_usable
    assert isinstance(result.rows[0], DirectObservationRowDiagnostic)
    assert result.rows[0].derivative_status == ("EXACT_PROVIDER_FIRST_DERIVATIVE")


def test_direct_observations_fail_before_native_evaluation_outside_domain() -> None:
    solvation_model = _fixed_water_factor_model(FIGIEL_WATER_SOLVATION_FACTOR_V1)
    solvation = _solvation_factor_problem(solvation_model)
    first = solvation.observations[0]
    wrong_pressure = _replace_observations(
        solvation,
        (
            replace(first, pressure_pa=200_000.0),
            *solvation.observations[1:],
        ),
    )
    low_molality = _replace_observations(
        solvation,
        (
            replace(
                first,
                formula_unit_molality_mol_per_kg=0.0009,
            ),
            *solvation.observations[1:],
        ),
    )
    with pytest.raises(ValueError, match="pressure or molality"):
        fit_parameters(wrong_pressure, solvation_model)
    with pytest.raises(ValueError, match="pressure or molality"):
        fit_parameters(low_molality, solvation_model)

    target = FIGIEL_BORN_DIAMETER_TRACER_V1.targets[0]
    born_model = _aqueous_model(target.component_order)
    born = _born_diameter_problem(born_model, 0)
    born_wrong_pressure = _replace_observations(
        born,
        (replace(born.observations[0], pressure_pa=200_000.0),),
    )
    with pytest.raises(ValueError, match="pressure is outside"):
        fit_parameters(born_wrong_pressure, born_model)


def test_installed_provider_advertises_scalar_pure_parameter_contracts() -> None:
    model = _pure_model()
    capabilities = tuple(
        capability
        for capability in parameter_capabilities(model)
        if not isinstance(capability, UnsupportedParameterCapability)
    )

    assert tuple(capability.family for capability in capabilities) == (
        ParameterFamily.SEGMENT_COUNT,
        ParameterFamily.SEGMENT_DIAMETER,
        ParameterFamily.DISPERSION_ENERGY_OVER_K,
    )
    assert tuple(capability.capability_id for capability in capabilities) == (
        "neutral_pure_segment_count_v1",
        "neutral_pure_segment_diameter_v1",
        "neutral_pure_dispersion_energy_over_k_v1",
    )
    assert tuple(capability.coordinate_kinds for capability in capabilities) == (
        ("amount", "volume", "segment_count"),
        ("amount", "volume", "segment_diameter"),
        ("amount", "volume", "dispersion_energy_over_k"),
    )
    assert tuple(capability.coordinate_units[-1] for capability in capabilities) == (
        "dimensionless",
        "angstrom",
        "kelvin",
    )
    assert all(capability.component_ids == ("methane",) for capability in capabilities)
    assert all(capability.identity_shape == "component" for capability in capabilities)
    assert all(
        capability.model_domain == "neutral_nonassociating_pure"
        for capability in capabilities
    )
    assert all(capability.state_coordinate_count == 2 for capability in capabilities)


@pytest.mark.campaign
def test_installed_provider_advertises_bounded_pure_association_contracts() -> None:
    capabilities = tuple(
        capability
        for capability in parameter_capabilities(_associating_pure_model())
        if not isinstance(capability, UnsupportedParameterCapability)
    )
    assert tuple(capability.family for capability in capabilities) == (
        ParameterFamily.ASSOCIATION_ENERGY_OVER_K,
        ParameterFamily.ASSOCIATION_VOLUME,
        ParameterFamily.PURE_ASSOCIATING_JOINT,
    )
    assert tuple(capability.capability_id for capability in capabilities) == (
        "neutral_pure_2b_association_energy_over_k_v1",
        "neutral_pure_2b_association_volume_v1",
        "neutral_pure_associating_joint_sigma_basis_v1",
    )
    assert tuple(capability.coordinate_kinds for capability in capabilities) == (
        ("amount", "volume", "association_energy_over_k"),
        ("amount", "volume", "association_volume"),
        (
            "amount",
            "volume",
            "segment_count",
            "segment_diameter",
            "dispersion_energy_over_k",
            "association_energy_over_k",
            "association_volume",
        ),
    )
    assert tuple(capability.coordinate_units[-1] for capability in capabilities) == (
        "kelvin",
        "dimensionless",
        "dimensionless",
    )
    assert all(capability.identity_shape == "model" for capability in capabilities)
    assert all(
        capability.model_domain == "neutral_associating_pure"
        for capability in capabilities
    )


def test_ordinary_sigma_2b_block_evaluates_exact_jacobian() -> None:
    sites = (("acceptor", 1), ("donor", 1))
    pairs = (("acceptor", "donor", 1500.0, 0.01),)
    model = _generic_associating_model(sites, pairs)
    capability = next(
        capability
        for capability in parameter_capabilities(model)
        if not isinstance(capability, UnsupportedParameterCapability)
        and capability.capability_id
        == "neutral_pure_associating_joint_sigma_basis_v1"
    )
    assert capability.active_parameter_count == 5
    assert capability.identity_shape == "model"
    problem = _generic_associating_problem(model, pairs)
    assert len(problem.parameters) == 5
    assert problem.parameter_slot_indices == tuple(range(5))
    assert all(
        isinstance(parameter.identity, ModelParameterIdentity)
        for parameter in problem.parameters[3:]
    )
    variables = (0.0,) * (len(problem.parameters) + 1)
    residuals, jacobian = _evaluate_parameters(problem, model, variables)
    assert len(residuals) == 2
    assert len(jacobian) == 2 * len(variables)
    assert all(math.isfinite(value) for value in (*residuals, *jacobian))
    step = 1.0e-6
    for column in range(len(variables)):
        lower = list(variables)
        upper = list(variables)
        lower[column] -= step
        upper[column] += step
        lower_residuals, _ = _evaluate_parameters(problem, model, tuple(lower))
        upper_residuals, _ = _evaluate_parameters(problem, model, tuple(upper))
        for row, (lower_value, upper_value) in enumerate(
            zip(lower_residuals, upper_residuals, strict=True)
        ):
            finite_difference = (upper_value - lower_value) / (2.0 * step)
            assert jacobian[row * len(variables) + column] == pytest.approx(
                finite_difference, rel=2.0e-5, abs=2.0e-6
            )

    source_row = problem.observations[0]
    parity_volume = 9.0e-5
    parity_density = 0.088 / parity_volume
    parity_rows = []
    for temperature in (300.0, 325.0, 350.0, 375.0, 400.0):
        public_state = model.state(
            T=temperature * unit_registry.kelvin,
            rho=(1.0 / parity_volume)
            * unit_registry.mole
            / unit_registry.meter**3,
            x=(1.0,),
        )
        public_pressure = float(public_state.pressure.to("pascal").magnitude)
        parity_rows.append(
            replace(
                source_row,
                row_id=f"public-state-parity-{temperature:g}-k",
                temperature_k=temperature,
                pressure_pa=public_pressure,
                pressure_scale_pa=public_pressure,
                density_kg_per_m3=parity_density,
                density_scale_kg_per_m3=parity_density,
                volume_origin_m3_per_mol=parity_volume,
                volume_start_m3_per_mol=parity_volume,
                volume_bounds_m3_per_mol=(8.0e-5, 1.0e-4),
            )
        )
    parity_observations = tuple(parity_rows)
    parity_problem = replace(
        problem,
        sources=(
            replace(
                problem.sources[0],
                canonical_dataset_sha256=canonical_dataset_sha256(
                    parity_observations
                ),
            ),
        ),
        observations=parity_observations,
        maximum_iterations=100,
    )
    parity_result = fit_parameters(parity_problem, model)
    assert parity_result.solver_converged
    assert parity_result.numerically_converged
    assert parity_result.workflow_valid
    assert (
        parity_result.jacobian.full_rank
        == parity_result.jacobian.variable_count
        == 10
    )
    assert parity_result.jacobian.projected_parameter_rank == 5
    fitted = tuple(parameter.final for parameter in parity_result.parameters)
    assert fitted == (
        3.2,
        3.5,
        280.0,
        1500.0,
        0.01,
    )
    assert max(
        abs(value)
        for row in parity_result.rows
        for value in row.scaled_residuals
    ) < 1.0e-12
    assert (
        parity_result.scientific_status
        == "NOT_ADJUDICATED_NO_APPROVED_SCIENTIFIC_CUTOFF"
    )
    replay_model = _generic_associating_model(
        sites,
        ((pairs[0][0], pairs[0][1], fitted[3], fitted[4]),),
        segment_count=fitted[0],
        segment_diameter_angstrom=fitted[1],
        dispersion_energy_over_k=fitted[2],
    )
    for row in parity_observations:
        replay_state = replay_model.state(
            T=row.temperature_k * unit_registry.kelvin,
            rho=(1.0 / parity_volume)
            * unit_registry.mole
            / unit_registry.meter**3,
            x=(1.0,),
        )
        replay_pressure = float(replay_state.pressure.to("pascal").magnitude)
        assert replay_pressure == pytest.approx(row.pressure_pa, rel=1.0e-13)


def test_generic_association_block_accepts_combined_saturation_rows() -> None:
    model = _generic_associating_model(
        (("acceptor", 1), ("donor", 1)),
        (("acceptor", "donor", 1500.0, 0.01),),
    )
    problem = _generic_associating_problem(
        model, (("acceptor", "donor", 1500.0, 0.01),)
    )
    row = PureSaturationObservation(
        row_id="manufactured-saturation",
        source_id=problem.sources[0].source_id,
        source_locator="manufactured combined saturation row",
        component_id="test-amine",
        temperature_k=350.0,
        pressure_pa=100_000.0,
        liquid_density_kg_per_m3=800.0,
        molar_mass_kg_per_mol=0.088,
        pressure_scale_pa=100_000.0,
        chemical_potential_scale=1.0,
        liquid_density_scale_kg_per_m3=800.0,
        liquid_volume_origin_m3_per_mol=0.088 / 800.0,
        liquid_volume_start_m3_per_mol=0.088 / 800.0,
        liquid_volume_bounds_m3_per_mol=(5.0e-5, 2.0e-4),
        vapor_volume_origin_m3_per_mol=0.03,
        vapor_volume_start_m3_per_mol=0.03,
        vapor_volume_bounds_m3_per_mol=(2.1e-4, 1.0),
        partition=ObservationPartition.TRAINING,
    )
    problem = replace(
        problem,
        sources=(
            replace(
                problem.sources[0],
                canonical_dataset_sha256=canonical_dataset_sha256((row,)),
            ),
        ),
        observations=(row,),
    )
    variables = (0.0,) * (len(problem.parameters) + 2)
    residuals, jacobian = _evaluate_parameters(problem, model, variables)

    assert len(residuals) == 4
    assert len(jacobian) == 4 * len(variables)
    assert all(math.isfinite(value) for value in (*residuals, *jacobian))
    step = 1.0e-6
    for column in range(len(variables)):
        lower = list(variables)
        upper = list(variables)
        lower[column] -= step
        upper[column] += step
        lower_residuals, _ = _evaluate_parameters(problem, model, tuple(lower))
        upper_residuals, _ = _evaluate_parameters(problem, model, tuple(upper))
        for row, (lower_value, upper_value) in enumerate(
            zip(lower_residuals, upper_residuals, strict=True)
        ):
            finite_difference = (upper_value - lower_value) / (2.0 * step)
            assert jacobian[row * len(variables) + column] == pytest.approx(
                finite_difference, rel=2.0e-5, abs=2.0e-6
            )


@pytest.mark.campaign
@pytest.mark.parametrize(
    "family",
    (
        ParameterFamily.ASSOCIATION_ENERGY_OVER_K,
        ParameterFamily.ASSOCIATION_VOLUME,
    ),
)
def test_pure_density_association_surface_reports_nonconverged_rank_diagnostics(
    family: ParameterFamily,
) -> None:
    model = _associating_pure_model()
    result = fit_parameters(_pure_density_problem(model, family), model)
    assert result.solution_usable
    assert result.jacobian.residual_count == 2
    assert result.jacobian.variable_count == 2
    assert result.jacobian.full_rank == 2
    assert result.jacobian.projected_parameter_rank == 1
    assert result.parameters[0].active_bound is None
    assert math.isfinite(result.parameters[0].final)
    assert not result.solver_converged
    assert not result.numerically_converged
    assert isinstance(result.rows[0], PureDensityRowDiagnostic)
    assert result.rows[0].evaluated
    assert all(math.isfinite(value) for value in result.rows[0].scaled_residuals)
    assert result.scientific_status == "NOT_ADJUDICATED_NO_APPROVED_SCIENTIFIC_CUTOFF"


@pytest.mark.parametrize(
    "family",
    (
        ParameterFamily.SEGMENT_COUNT,
        ParameterFamily.SEGMENT_DIAMETER,
        ParameterFamily.DISPERSION_ENERGY_OVER_K,
    ),
)
def test_general_engine_fits_one_pure_component_parameter(
    family: ParameterFamily,
) -> None:
    model = _pure_model()
    result = fit_parameters(_pure_problem(model, family, all_training_rows=True), model)

    assert result.solver_converged
    assert result.numerically_converged
    assert result.workflow_valid
    assert result.parameters[0].family is family
    assert result.parameters[0].active_bound is None
    assert result.jacobian.full_rank == result.jacobian.variable_count == 9
    assert result.jacobian.projected_parameter_rank == 1
    expected = {
        ParameterFamily.SEGMENT_COUNT: 1.0001569260577763,
        ParameterFamily.SEGMENT_DIAMETER: 3.7063548743836034,
        ParameterFamily.DISPERSION_ENERGY_OVER_K: 150.00325287725062,
    }
    assert result.parameters[0].final == pytest.approx(expected[family], abs=2.0e-10)
    assert all(
        isinstance(row, PureSaturationRowDiagnostic) and row.evaluated
        for row in result.rows
    )


def test_general_engine_fits_joint_pure_parameter_vector() -> None:
    model = _pure_model()
    result = fit_parameters(_joint_pure_problem(model), model)

    assert result.solver_converged
    assert result.numerically_converged
    assert result.workflow_valid
    assert tuple(parameter.family for parameter in result.parameters) == (
        ParameterFamily.SEGMENT_COUNT,
        ParameterFamily.SEGMENT_DIAMETER,
        ParameterFamily.DISPERSION_ENERGY_OVER_K,
    )
    assert tuple(parameter.final for parameter in result.parameters) == pytest.approx(
        (0.9932081279826167, 3.717121437945618, 150.4888402511307),
        rel=2.0e-9,
        abs=2.0e-9,
    )
    assert not any(parameter.active_bound for parameter in result.parameters)
    assert result.jacobian.residual_count == 16
    assert result.jacobian.variable_count == 11
    assert result.jacobian.full_rank == 11
    assert result.jacobian.projected_parameter_rank == 3
    assert result.confirmation_count == 1
    assert result.confirmations_usable
    assert all(
        isinstance(row, PureSaturationRowDiagnostic) and row.evaluated
        for row in result.rows
    )


def test_public_preparation_preserves_joint_pure_problem_semantics() -> None:
    model = _pure_model()
    direct = _joint_pure_problem(model)

    prepared = _prepare_existing_problem(model, direct, "pure_saturation")

    assert prepared.problem == direct
    assert prepared.preflight().ready


def test_public_preparation_preserves_fixed_2b_problem_semantics() -> None:
    pairs = (("acceptor", "donor", 1500.0, 0.01),)
    model = _generic_associating_model(
        (("acceptor", 1), ("donor", 1)),
        pairs,
    )
    direct = _generic_associating_problem(model, pairs)

    prepared = _prepare_existing_problem(model, direct, "pure_density")

    assert prepared.problem == direct
    assert tuple(
        coordinate.family for coordinate in prepared.problem.parameters
    ) == (
        ParameterFamily.SEGMENT_COUNT,
        ParameterFamily.SEGMENT_DIAMETER,
        ParameterFamily.DISPERSION_ENERGY_OVER_K,
        ParameterFamily.ASSOCIATION_ENERGY_OVER_K,
        ParameterFamily.ASSOCIATION_VOLUME,
    )


def test_public_preparation_preserves_direct_observable_problem_semantics() -> None:
    target = FIGIEL_BORN_DIAMETER_TRACER_V1.targets[0]
    model = _aqueous_model(target.component_order)
    direct = _born_diameter_problem(model, 0)

    prepared = _prepare_existing_problem(model, direct, "solvation_gibbs")

    assert prepared.problem == direct
    assert prepared.preflight().ready


def test_native_joint_pure_adapter_rejects_reordered_slots() -> None:
    model = _pure_model()
    problem = _joint_pure_problem(model)
    capabilities = tuple(
        _capability(model, parameter.family) for parameter in problem.parameters
    )
    payload = list(_native_payload(problem, capabilities[0]))
    payload[-1] = (0, 2, 1)

    with pytest.raises(RuntimeError, match="declared m, sigma"):
        parameter_regression._native.evaluate_general(
            native_sdk(model), tuple(payload), (0.0,) * 11
        )


def test_native_general_engine_rejects_nonpositive_pure_residual_scale() -> None:
    model = _pure_model()
    problem = _pure_problem(model, ParameterFamily.SEGMENT_COUNT)
    capability = _capability(model, ParameterFamily.SEGMENT_COUNT)
    payload = list(_native_payload(problem, capability))
    rows = list(payload[-3])
    row = list(rows[0])
    row[4] = 0.0
    rows[0] = tuple(row)
    payload[-3] = tuple(rows)

    with pytest.raises(RuntimeError, match="positive and ordered"):
        parameter_regression._native.evaluate_general(
            native_sdk(model), tuple(payload), (0.0, 0.0, 0.0)
        )


def test_native_general_engine_rejects_nonfinite_scaled_result() -> None:
    model = _pure_model()
    problem = _pure_problem(model, ParameterFamily.SEGMENT_COUNT)
    capability = _capability(model, ParameterFamily.SEGMENT_COUNT)
    payload = list(_native_payload(problem, capability))
    rows = list(payload[-3])
    row = list(rows[0])
    row[4] = 1.0e-320
    rows[0] = tuple(row)
    payload[-3] = tuple(rows)

    with pytest.raises(RuntimeError, match="assembled residual or Jacobian"):
        parameter_regression._native.evaluate_general(
            native_sdk(model), tuple(payload), (0.0, 0.0, 0.0)
        )


def test_general_kij_fit_reports_rank_confirmation_and_partition_isolation() -> None:
    model = _model()
    training = _row()
    held = replace(
        _row("may-2015-held", ObservationPartition.HELD_OUT),
        pressure_pa=2.2e6,
        pressure_scale_pa=2.2e6,
    )

    problem = _problem(model, (training, held))
    result = fit_parameters(problem, model)
    perturbed_held = replace(held, pressure_pa=8.0e6, pressure_scale_pa=8.0e6)
    repeated = fit_parameters(_problem(model, (training, perturbed_held)), model)

    assert isinstance(result, RegressionResult)
    assert result.problem == problem
    assert result.capabilities == (parameter_capabilities(model)[0],)
    assert result.parameters[0].final == repeated.parameters[0].final
    assert result.final_cost == repeated.final_cost
    assert (
        result.parameters[0].lower_bound
        <= result.parameters[0].final
        <= result.parameters[0].upper_bound
    )
    assert result.parameters[0].transform_origin == 0.0
    assert result.parameters[0].transform_scale == 0.01
    assert result.jacobian.full_rank == 3
    assert result.jacobian.projected_parameter_rank == 1
    assert result.confirmation_count == 2
    assert result.training_row_count == 1
    assert result.held_out_row_count == 1
    assert result.stress_row_count == 0
    assert tuple(row.partition for row in result.rows) == ("training", "held_out")
    assert result.rows[0].observed_pressure_pa == training.pressure_pa
    assert result.rows[0].liquid_model_pressure_pa == pytest.approx(
        training.pressure_pa
        + result.rows[0].scaled_residuals[0] * training.pressure_scale_pa
    )
    assert result.rows[0].derivative_status == "EXACT_PROVIDER_HESSIAN"
    assert result.rows[0].status == "evaluated"
    assert result.evaluated_row_count == 2
    assert result.skipped_row_count == 0
    assert result.failed_row_count == 0
    assert result.predictive_status == "NOT_ADJUDICATED_NO_APPROVED_HELD_OUT_CUTOFF"
    assert result.residual_evaluation_count > 0
    assert type(result.residual_evaluation_count) is int
    assert result.jacobian_evaluation_count > 0
    assert type(result.jacobian_evaluation_count) is int


def test_general_engine_tiny_solver_budget_stops_without_numerical_convergence() -> (
    None
):
    model = _model()
    problem = replace(_problem(model), maximum_solver_time_seconds=1.0e-12)

    result = fit_parameters(problem, model)

    assert result.termination != "CONVERGENCE"
    assert not result.numerically_converged


def test_general_lij_fit_reuses_the_exact_lifted_pair_engine() -> None:
    model = _model()

    result = fit_parameters(
        _problem(model, family=ParameterFamily.L_IJ),
        model,
    )

    assert result.parameters[0].family is ParameterFamily.L_IJ
    assert result.capabilities[0].family is ParameterFamily.L_IJ
    assert result.jacobian.residual_count == 4
    assert result.jacobian.variable_count == 3
    assert result.jacobian.projected_parameter_rank == 1
    assert result.parameters[0].active_bound is None
    assert result.confirmation_count == 2
    assert result.rows[0].derivative_status == "EXACT_PROVIDER_HESSIAN"


@pytest.mark.parametrize(
    ("diagnostic_index", "diagnostic_value"),
    (
        pytest.param(11, 2, id="full-rank-deficient"),
        pytest.param(14, 0, id="projected-rank-zero"),
        pytest.param(12, 1.0e11, id="full-conditioning-over-limit"),
        pytest.param(15, 1.0e11, id="projected-conditioning-over-limit"),
    ),
)
def test_numerical_status_requires_complete_rank_and_conditioning(
    monkeypatch: pytest.MonkeyPatch,
    diagnostic_index: int,
    diagnostic_value: float,
) -> None:
    model = _model()
    problem = _problem(model)
    native = _mock_general_native_result()
    monkeypatch.setattr(
        parameter_regression._native, "solve_general", lambda *_: native
    )

    converged = fit_parameters(problem, model)
    assert converged.solver_converged
    assert converged.numerically_converged
    assert converged.jacobian.full_rank == 3

    degraded = list(native)
    degraded[diagnostic_index] = diagnostic_value
    monkeypatch.setattr(
        parameter_regression._native,
        "solve_general",
        lambda *_: tuple(degraded),
    )
    diagnostically_stopped = fit_parameters(problem, model)
    assert diagnostically_stopped.solver_converged
    assert not diagnostically_stopped.numerically_converged

    monkeypatch.setattr(
        parameter_regression._native,
        "solve_general",
        lambda *_: ("NO_CONVERGENCE", *native[1:]),
    )
    stopped = fit_parameters(problem, model)
    assert not stopped.solver_converged
    assert not stopped.numerically_converged


def test_active_bound_diagnostic_is_preserved(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    model = _model()
    problem = _problem(model)
    native = list(_mock_general_native_result())
    native[5] = problem.parameters[0].upper_bound
    native[6] = 0.0
    native[7] = "upper"
    monkeypatch.setattr(
        parameter_regression._native,
        "solve_general",
        lambda *_: tuple(native),
    )

    result = fit_parameters(problem, model)

    assert result.parameters[0].final == problem.parameters[0].upper_bound
    assert result.parameters[0].active_bound == "upper"
    assert result.parameters[0].active_bound_distance == 0.0


def test_rows_outside_provider_temperature_domain_fail_before_ceres(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    model = _model()
    capability = _capability(model)
    outside = replace(_row(), temperature_k=capability.temperature_max_k + 1.0)
    problem = _problem(model, (outside,))
    monkeypatch.setattr(
        parameter_regression._native,
        "solve_general",
        lambda *_: pytest.fail("Ceres must not start outside the Provider domain"),
    )

    with pytest.raises(ValueError, match="temperature"):
        fit_parameters(problem, model)


def test_parameter_fingerprint_mismatch_fails_before_ceres(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    model = _model()
    problem = _problem(model)
    mismatched = replace(
        problem,
        parameters=(
            replace(
                problem.parameters[0],
                provider_parameter_fingerprint=f"sha256:{'0' * 64}",
            ),
        ),
    )
    monkeypatch.setattr(
        parameter_regression._native,
        "solve_general",
        lambda *_: pytest.fail("Ceres must not start for a fingerprint mismatch"),
    )

    with pytest.raises(ValueError, match="does not match"):
        fit_parameters(mismatched, model)


def test_phase_observation_contract_mismatch_fails_before_ceres(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    model = _model()
    problem = _problem(model)
    capability = replace(
        _capability(model),
        observation_contract="aqueous_mean_ionic_activity",
    )
    monkeypatch.setattr(
        parameter_regression,
        "parameter_capabilities",
        lambda _: (capability,),
    )
    monkeypatch.setattr(
        parameter_regression._native,
        "solve_general",
        lambda *_: pytest.fail(
            "Ceres must not start for an observation-contract mismatch"
        ),
    )

    with pytest.raises(ValueError, match="phase observation"):
        fit_parameters(problem, model)


def test_provider_failure_returns_diagnostic_result() -> None:
    model = _model()
    invalid = replace(
        _row(),
        liquid_volume_origin_m3_per_mol=1.0e-10,
        liquid_volume_start_m3_per_mol=1.0e-10,
        liquid_volume_bounds_m3_per_mol=(5.0e-11, 2.0e-10),
    )

    result = fit_parameters(_problem(model, (invalid,)), model)

    assert not result.solver_converged
    assert not result.numerically_converged
    assert not result.workflow_valid
    assert result.failure_reasons
    assert not result.rows[0].evaluated
    assert result.rows[0].derivative_status == "UNAVAILABLE"
    assert result.rows[0].status == "failed"
    assert result.evaluated_row_count == 0
    assert result.failed_row_count == 1


def _audited_may_rows() -> tuple[FixedCompositionVleObservation, ...]:
    data_path = (
        Path(__file__).parents[1] / "evidence" / "may-2015-methane-ethane-vle.csv"
    )
    with data_path.open(newline="", encoding="utf-8") as stream:
        records = tuple(csv.DictReader(stream))
    gas_constant = 8.31446261815324
    return tuple(
        FixedCompositionVleObservation(
            row_id=record["row_id"],
            source_id="may-2015",
            source_locator=f"{data_path.name}:{record['row_id']}",
            component_ids=("methane", "ethane"),
            temperature_k=float(record["T_K"]),
            pressure_pa=float(record["P_Pa"]),
            liquid_mole_fraction_first=float(record["x_methane"]),
            vapor_mole_fraction_first=float(record["y_methane"]),
            pressure_scale_pa=float(record["P_Pa"]),
            chemical_potential_scales=(1.0, 1.0),
            liquid_volume_origin_m3_per_mol=6.5e-5,
            liquid_volume_start_m3_per_mol=6.5e-5,
            liquid_volume_bounds_m3_per_mol=(2.0e-5, 1.0e-4),
            vapor_volume_origin_m3_per_mol=(
                gas_constant * float(record["T_K"]) / float(record["P_Pa"])
            ),
            vapor_volume_start_m3_per_mol=(
                gas_constant * float(record["T_K"]) / float(record["P_Pa"])
            ),
            vapor_volume_bounds_m3_per_mol=(1.0e-4, 1.0e-2),
            partition=ObservationPartition.TRAINING,
        )
        for record in records
    )


def _may_methane_propane_records() -> tuple[dict[str, str], ...]:
    data_path = (
        Path(__file__).parents[1] / "evidence" / "may-2015-methane-propane-vle.csv"
    )
    with data_path.open(newline="", encoding="utf-8") as stream:
        return tuple(csv.DictReader(stream))


def _may_methane_propane_rows() -> tuple[FixedCompositionVleObservation, ...]:
    gas_constant = 8.31446261815324
    data_path = (
        Path(__file__).parents[1] / "evidence" / "may-2015-methane-propane-vle.csv"
    )
    return tuple(
        FixedCompositionVleObservation(
            row_id=record["row_id"],
            source_id="may-2015-methane-propane",
            source_locator=f"{data_path.name}:{record['row_id']}",
            component_ids=("methane", "propane"),
            temperature_k=float(record["T_K"]),
            pressure_pa=1000.0 * float(record["p_kPa"]),
            liquid_mole_fraction_first=float(record["x_methane"]),
            vapor_mole_fraction_first=1.0 - float(record["y_propane"]),
            pressure_scale_pa=1000.0 * float(record["p_kPa"]),
            chemical_potential_scales=(1.0, 1.0),
            liquid_volume_origin_m3_per_mol=4.0e-5,
            liquid_volume_start_m3_per_mol=4.0e-5,
            liquid_volume_bounds_m3_per_mol=(3.0e-5, 2.0e-4),
            vapor_volume_origin_m3_per_mol=(
                gas_constant * float(record["T_K"]) / (1000.0 * float(record["p_kPa"]))
            ),
            vapor_volume_start_m3_per_mol=(
                gas_constant * float(record["T_K"]) / (1000.0 * float(record["p_kPa"]))
            ),
            vapor_volume_bounds_m3_per_mol=(5.0e-5, 1.0e-2),
            partition=ObservationPartition.TRAINING,
        )
        for record in _may_methane_propane_records()
    )


def _may_methane_propane_problem(
    model: Mixture,
    rows: tuple[FixedCompositionVleObservation, ...] | None = None,
) -> RegressionProblem:
    observations = rows or _may_methane_propane_rows()
    capability = next(
        (
            candidate
            for candidate in parameter_capabilities(model)
            if not isinstance(candidate, UnsupportedParameterCapability)
            and candidate.capability_id == "neutral_binary_phase_kij_v1"
        ),
        None,
    )
    assert capability is not None
    assert capability.family is ParameterFamily.K_IJ
    assert capability.observation_contract == "fixed_composition_helmholtz_phase"
    assert capability.component_ids == ("methane", "propane")
    assert capability.active_component_ids == ("methane", "propane")
    assert capability.model_domain == "neutral_nonassociating_binary"
    assert all(row.component_ids == capability.component_ids for row in observations)
    source = SourceDescriptor(
        source_id="may-2015-methane-propane",
        citation=(
            "May et al. (2015), Reference Quality Vapor--Liquid Equilibrium "
            "Data for the Binary Systems Methane + Ethane, + Propane, "
            "J. Chem. Eng. Data 60, 3606--3620, DOI 10.1021/acs.jced.5b00610."
        ),
        durable_locator="evidence/may-2015-methane-propane-vle.csv",
        source_artifact_sha256=(
            "53fd1bdd55dc6807ec76cf88626438d8dfceb3ec09149d4405ea36cfbe6b842a"
        ),
        canonical_dataset_sha256=canonical_dataset_sha256(observations),
        transformation_record=(
            "Transcribed Table 6 values unchanged; converted p/kPa to P/Pa "
            "by multiplication by 1000 and derived methane vapor fraction as "
            "1 - y_propane. Source uncertainties remain descriptive inputs "
            "and are not fit weights or cutoffs."
        ),
        units_and_bases=(
            "Source T/K, p/kPa, mole fractions, and standard/combined "
            "uncertainties; model pressure/Pa and dimensionless mu/RT scales."
        ),
        use_basis=(
            "Direct experimental binary VLE observations used for one "
            "in-sample constant-k_ij reference fit."
        ),
        residual_scale_rationale=(
            "Pressure residuals use the observed pressure magnitude; both "
            "chemical-potential residuals use unit mu/RT scales."
        ),
    )
    parameter = ParameterCoordinate(
        family=ParameterFamily.K_IJ,
        identity=PairParameterIdentity("methane", "propane"),
        capability_id=capability.capability_id,
        provider_parameter_fingerprint=capability.parameter_fingerprint,
        provider_topology_fingerprint=capability.topology_fingerprint,
        unit="1",
        transform=AffineParameterTransform(origin=0.0, scale=0.01),
        lower_bound=-0.15,
        upper_bound=0.10,
    )
    return RegressionProblem(
        sources=(source,),
        parameters=(parameter,),
        parameter_slot_indices=(0,),
        start_vectors=((0.0,), (-0.05,), (0.05,)),
        observations=observations,
        maximum_condition_number=1.0e10,
        maximum_iterations=100,
        maximum_solver_time_seconds=180.0,
        function_tolerance=1.0e-12,
        gradient_tolerance=1.0e-12,
        parameter_tolerance=1.0e-12,
        confirmation_parameter_scaled_max_delta=1.0e-5,
        confirmation_cost_relative_delta=1.0e-8,
    )


def test_may_methane_propane_source_identity_and_transform() -> None:
    data_path = (
        Path(__file__).parents[1] / "evidence" / "may-2015-methane-propane-vle.csv"
    )
    assert (
        canonical_dataset_sha256(_may_methane_propane_rows())
        == _MAY_METHANE_PROPANE_CANONICAL_SHA256
    )
    assert (
        data_path.read_bytes()
        and __import__("hashlib").sha256(data_path.read_bytes()).hexdigest()
        == "97a07b274dc4da6a281614f3fd39c520ebd6678776413746b13bc8665113c529"
    )
    records = _may_methane_propane_records()
    assert len(records) == 22
    assert len({record["row_id"] for record in records}) == 22
    expected = (
        (283.38, 7630, 0.4586, 0.0010, 0.0018, 0.2045, 0.0006, 0.0009),
        (273.48, 7102, 0.4589, 0.0010, 0.0018, 0.1614, 0.0003, 0.0007),
        (263.51, 6597, 0.4623, 0.0007, 0.0018, 0.1197, 0.0012, 0.0013),
        (263.56, 6560, 0.4601, 0.0013, 0.0021, 0.1174, 0.0018, 0.0019),
        (253.59, 6100, 0.4708, 0.0011, 0.0021, 0.0947, 0.0004, 0.0006),
        (243.61, 5611, 0.4857, 0.0015, 0.0024, 0.0715, 0.0001, 0.0004),
        (243.64, 5635, 0.4907, 0.0010, 0.0022, 0.0719, 0.0001, 0.0004),
        (243.62, 5544, 0.4832, 0.0014, 0.0024, 0.0715, 0.0002, 0.0004),
        (233.67, 5119, 0.5077, 0.0012, 0.0026, 0.0518, 0.0003, 0.0004),
        (233.61, 5077, 0.5037, 0.0021, 0.0031, 0.0520, 0.0003, 0.0004),
        (223.54, 4622, 0.5333, 0.0006, 0.0027, 0.0368, 0.0005, 0.0005),
        (223.36, 4628, 0.5331, 0.0012, 0.0030, 0.0366, 0.0002, 0.0003),
        (213.38, 4108, 0.5614, 0.0005, 0.0033, 0.0245, 0.0003, 0.0003),
        (203.40, 3576, 0.5952, 0.0014, 0.0044, 0.0162, 0.0005, 0.0005),
        (243.62, 891, 0.0710, 0.0005, 0.0019, 0.2076, 0.0003, 0.0041),
        (243.62, 3909, 0.3442, 0.0003, 0.0019, 0.0743, 0.0004, 0.0006),
        (243.63, 5500, 0.4788, 0.0014, 0.0024, 0.0704, 0.0002, 0.0004),
        (243.61, 5544, 0.4764, 0.0015, 0.0025, 0.0734, 0.0002, 0.0004),
        (243.62, 6402, 0.5565, 0.0006, 0.0022, 0.0726, 0.0006, 0.0007),
        (243.63, 6989, 0.6082, 0.0015, 0.0026, 0.0768, 0.0005, 0.0006),
        (243.62, 7530, 0.6572, 0.0010, 0.0025, 0.0779, 0.0016, 0.0016),
        (243.62, 7943, 0.6961, 0.0013, 0.0027, 0.0855, 0.0005, 0.0007),
    )
    observed = tuple(
        (
            float(record["T_K"]),
            int(record["p_kPa"]),
            float(record["x_methane"]),
            float(record["u_x_methane"]),
            float(record["uc_x_methane"]),
            float(record["y_propane"]),
            float(record["u_y_propane"]),
            float(record["uc_y_propane"]),
        )
        for record in records
    )
    assert observed == expected
    rows = _may_methane_propane_rows()
    problem = _may_methane_propane_problem(_methane_propane_model(), rows)
    parameters = _methane_propane_parameters()
    pair_records = tuple(
        record
        for record in parameters.records
        if isinstance(record, PairParameterRecord)
    )
    assert parameters.bundle_purpose == "user-provided"
    assert problem.observations[0].component_ids == ("methane", "propane")
    assert len(pair_records) == 1
    assert pair_records[0].component_id_a == "methane"
    assert pair_records[0].component_id_b == "propane"
    assert float(pair_records[0].value) == 0.0
    assert problem.sources[0].source_artifact_sha256 == (
        "53fd1bdd55dc6807ec76cf88626438d8dfceb3ec09149d4405ea36cfbe6b842a"
    )
    assert tuple(row.pressure_pa for row in rows) == tuple(
        1000.0 * values[1] for values in expected
    )
    assert tuple(row.vapor_mole_fraction_first for row in rows) == tuple(
        1.0 - values[5] for values in expected
    )
    assert all(
        record["u_x_methane"]
        and record["uc_x_methane"]
        and record["u_y_propane"]
        and record["uc_y_propane"]
        for record in records
    )


def test_may_methane_propane_kij_is_invariant_to_row_order_and_pair_identity() -> None:
    # Two source rows keep the invariance check independent of the retained
    # 22-row campaign while still exercising a nontrivial row reversal.
    rows = _may_methane_propane_rows()[:2]
    base_model = _methane_propane_model()
    base_result = fit_parameters(
        _may_methane_propane_problem(base_model, rows), base_model
    )

    reversed_model = _methane_propane_model()
    reversed_result = fit_parameters(
        _may_methane_propane_problem(reversed_model, tuple(reversed(rows))),
        reversed_model,
    )

    identity_model = _methane_propane_model()
    identity_problem = _may_methane_propane_problem(identity_model, rows)
    identity_problem = replace(
        identity_problem,
        parameters=(
            replace(
                identity_problem.parameters[0],
                identity=PairParameterIdentity("propane", "methane"),
            ),
        ),
    )
    identity_result = fit_parameters(identity_problem, identity_model)

    base_signature = _general_result_signature(base_result)
    variants = (reversed_result, identity_result)
    assert base_result.solver_converged
    assert base_result.numerically_converged
    assert base_result.workflow_valid
    for variant in variants:
        assert variant.solver_converged
        assert variant.numerically_converged
        assert variant.workflow_valid
        signature = _general_result_signature(variant)
        assert signature[0] == pytest.approx(base_signature[0], rel=0.0, abs=1.0e-14)
        assert signature[1] == pytest.approx(base_signature[1], rel=0.0, abs=1.0e-14)
        assert signature[2:] == base_signature[2:]

    assert identity_problem.parameters[0].identity.component_ids == (
        "methane",
        "propane",
    )
    assert identity_result.parameters[0].component_ids == (
        "methane",
        "propane",
    )


@pytest.mark.campaign
def test_may_methane_propane_campaign_reproduces_reference_fit() -> None:
    model = _methane_propane_model()
    problem = _may_methane_propane_problem(model)

    result = fit_parameters(problem, model)
    repeat = fit_parameters(problem, _methane_propane_model())

    assert len(problem.observations) == 22
    assert problem.parameters[0].identity == PairParameterIdentity("methane", "propane")
    assert problem.start_vectors == ((0.0,), (-0.05,), (0.05,))
    assert result.solver_converged
    assert result.numerically_converged
    assert result.workflow_valid
    assert result.solution_usable
    assert result.termination == "CONVERGENCE"
    assert result.jacobian.residual_count == 88
    assert result.jacobian.variable_count == 45
    assert result.jacobian.full_rank == 45
    assert result.jacobian.projected_parameter_rank == 1
    assert result.parameters[0].active_bound is None
    assert result.confirmation_count == 2
    assert result.confirmations_usable
    assert result.training_row_count == 22
    assert result.held_out_row_count == 0
    assert result.stress_row_count == 0
    assert result.evaluated_row_count == 22
    assert result.skipped_row_count == 0
    assert result.failed_row_count == 0
    assert result.physical_status == "NOT_ADJUDICATED_NO_ROW_ACCEPTANCE_CRITERIA"
    assert result.scientific_status == ("NOT_ADJUDICATED_NO_APPROVED_SCIENTIFIC_CUTOFF")
    assert result.predictive_status == "NOT_ADJUDICATED_NO_APPROVED_HELD_OUT_CUTOFF"
    assert result.residual_evaluation_count > 0
    assert isinstance(result.residual_evaluation_count, int)
    assert result.jacobian_evaluation_count > 0
    assert isinstance(result.jacobian_evaluation_count, int)
    assert all(
        row.derivative_status == "EXACT_PROVIDER_HESSIAN"
        and row.evaluated
        and row.status == "evaluated"
        and not row.failure_reason
        for row in result.rows
    )
    assert result.parameters[0].final == pytest.approx(
        0.0038919335722629794, rel=0.0, abs=1.0e-14
    )
    assert result.final_cost == pytest.approx(0.03734758119771876, rel=0.0, abs=1.0e-14)

    # A second complete run establishes the observed deterministic
    # repeatability; no scientific or wall-time cutoff is inferred here.
    repeatability_parameter_delta = abs(
        result.parameters[0].final - repeat.parameters[0].final
    )
    repeatability_cost_delta = abs(result.final_cost - repeat.final_cost)
    assert repeatability_parameter_delta <= 1.0e-14
    assert repeatability_cost_delta <= 1.0e-14
    assert result.parameters[0].final == repeat.parameters[0].final
    assert result.final_cost == repeat.final_cost


@pytest.mark.campaign
def test_all_audited_may_rows_reproduce_the_general_kij_reference_fit() -> None:
    model = _model()
    rows = _audited_may_rows()
    problem = _problem(model, rows)

    result = fit_parameters(problem, model)

    assert len(rows) == 17
    assert result.termination == "CONVERGENCE"
    assert result.solution_usable
    assert result.jacobian.residual_count == 68
    assert result.jacobian.variable_count == 35
    assert result.jacobian.full_rank == 35
    assert result.jacobian.projected_parameter_rank == 1
    assert result.parameters[0].active_bound is None
    assert result.parameters[0].final == pytest.approx(
        -0.00843032298906253, rel=0.0, abs=2.0e-12
    )
    assert result.confirmation_count == 2
    assert result.training_row_count == 17
    assert result.predictive_status == "NOT_ADJUDICATED_NO_APPROVED_HELD_OUT_CUTOFF"


@pytest.mark.campaign
def test_all_audited_may_rows_are_fit_ready_for_general_lij() -> None:
    model = _model()
    rows = _audited_may_rows()

    result = fit_parameters(
        _problem(model, rows, family=ParameterFamily.L_IJ),
        model,
    )

    assert result.termination == "CONVERGENCE"
    assert result.solution_usable
    assert result.jacobian.residual_count == 68
    assert result.jacobian.variable_count == 35
    assert result.jacobian.full_rank == 35
    assert result.jacobian.projected_parameter_rank == 1
    assert result.parameters[0].active_bound is None
    assert result.confirmation_count == 2
    assert result.training_row_count == 17
    assert result.failed_row_count == 0
    assert result.predictive_status == "NOT_ADJUDICATED_NO_APPROVED_HELD_OUT_CUTOFF"
