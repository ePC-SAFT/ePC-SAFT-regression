from __future__ import annotations

from dataclasses import dataclass
import csv
import hashlib
from importlib.resources import files
import io
import math


EXPECTED_HEADER = ("species", "T_K", "p_sat_Pa", "rho_sat_liq_kg_m3", "source")
SOURCE_RETRIEVED_ON = "2026-07-17"
SOURCE_USE_BASIS = (
    "NIST Standard Reference Data retained as compact source-backed candidate evidence; "
    "redistribution and use remain subject to the NIST SRD terms"
)
SOURCE_UNITS = (
    ("temperature", "K"),
    ("pressure", "Pa"),
    ("saturated_liquid_mass_density", "kg/m3"),
)

METHANE_DATA_SHA256 = "a5e16df3bf8ec78483fc340782cddc89ab8b284a9f6dfaecd6cda3ffde579227"
METHANE_PACKAGED_DATA_SHA256 = (
    "dec64d5a6cac414a4a92393a0d728fa27c02135c6a159d0d1881d7b6dde6d26c"
)
METHANE_SOURCE_ID = "nist-webbook-srd69-methane-saturation"
METHANE_SOURCE_URL = (
    "https://webbook.nist.gov/cgi/fluid.cgi?Action=Data&Wide=on&ID=C74828&"
    "Type=SatP&Digits=8&THigh=180&TLow=100&TInc=10&RefState=DEF&TUnit=K&"
    "PUnit=Pa&DUnit=kg%2Fm3&HUnit=kJ%2Fmol&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm"
)
METHANE_SOURCE_CITATION = (
    "NIST Chemistry WebBook, SRD 69, methane (CAS 74-82-8) fluid properties"
)
METHANE_SOURCE_LOCATOR = "Saturation properties query, 100 K through 180 K in 10 K increments"
METHANE_SOURCE_TRANSFORMATION = (
    "Exact retained CSV fields and decimal strings; CRLF line endings normalized to LF. "
    f"Retained source SHA-256: {METHANE_DATA_SHA256}; packaged SHA-256: "
    f"{METHANE_PACKAGED_DATA_SHA256}."
)
METHANE_TEMPERATURES_K = tuple(float(value) for value in range(100, 181, 10))
METHANE_TRAINING_TEMPERATURES_K = (110.0, 130.0, 150.0, 170.0)
METHANE_HELD_OUT_TEMPERATURES_K = (100.0, 120.0, 140.0, 160.0, 180.0)

ETHANE_DATA_SHA256 = "ed09b8781acfb7025ca505878b884f6353ddd9f3f4bd7aae2e6df88bbe847a67"
ETHANE_PACKAGED_DATA_SHA256 = (
    "b01333e827933c0a7148672c8ae3eef78393320c0d18f2c4d5a0fc40d9bef6b2"
)
ETHANE_SOURCE_ID = "nist-webbook-srd69-ethane-saturation"
ETHANE_SOURCE_URL = (
    "https://webbook.nist.gov/cgi/fluid.cgi?Action=Data&Wide=on&ID=C74840&"
    "Type=SatP&Digits=8&THigh=280&TLow=100&TInc=20&RefState=DEF&TUnit=K&"
    "PUnit=Pa&DUnit=kg%2Fm3&HUnit=kJ%2Fmol&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm"
)
ETHANE_SOURCE_CITATION = (
    "NIST Chemistry WebBook, SRD 69, ethane (CAS 74-84-0) fluid properties"
)
ETHANE_SOURCE_LOCATOR = "Saturation properties query, 100 K through 280 K in 20 K increments"
ETHANE_SOURCE_TRANSFORMATION = (
    "Exact retained CSV fields and decimal strings; CRLF line endings normalized to LF. "
    f"Retained source SHA-256: {ETHANE_DATA_SHA256}; packaged SHA-256: "
    f"{ETHANE_PACKAGED_DATA_SHA256}."
)
ETHANE_TEMPERATURES_K = tuple(float(value) for value in range(100, 281, 20))
ETHANE_TRAINING_TEMPERATURES_K = (140.0, 180.0, 220.0, 260.0)
ETHANE_HELD_OUT_TEMPERATURES_K = (120.0, 160.0, 200.0, 240.0)
ETHANE_STRESS_TEMPERATURES_K = (100.0, 280.0)

# Retained import names used only by the accepted methane receipt tooling.
EXPECTED_DATA_SHA256 = METHANE_DATA_SHA256
EXPECTED_PACKAGED_DATA_SHA256 = METHANE_PACKAGED_DATA_SHA256


def _positive_finite(value: float, field: str) -> None:
    if not math.isfinite(value):
        raise ValueError(f"{field} must be finite")
    if value <= 0.0:
        raise ValueError(f"{field} must be positive")


def _source_fields(source_id: str) -> tuple[str, str, str, str, str, str]:
    if source_id == METHANE_SOURCE_ID:
        return (
            METHANE_SOURCE_CITATION,
            METHANE_SOURCE_LOCATOR,
            METHANE_SOURCE_URL,
            METHANE_SOURCE_TRANSFORMATION,
            METHANE_DATA_SHA256,
            METHANE_PACKAGED_DATA_SHA256,
        )
    if source_id == ETHANE_SOURCE_ID:
        return (
            ETHANE_SOURCE_CITATION,
            ETHANE_SOURCE_LOCATOR,
            ETHANE_SOURCE_URL,
            ETHANE_SOURCE_TRANSFORMATION,
            ETHANE_DATA_SHA256,
            ETHANE_PACKAGED_DATA_SHA256,
        )
    raise ValueError("source_id must identify an admitted pure-saturation table")


@dataclass(frozen=True, slots=True)
class SourceIdentity:
    source_id: str
    citation: str
    locator: str
    url: str
    retrieved_on: str
    use_basis: str
    transformation: str
    units: tuple[tuple[str, str], ...]
    data_sha256: str
    packaged_data_sha256: str

    def __post_init__(self) -> None:
        for field in (
            "source_id",
            "citation",
            "locator",
            "url",
            "retrieved_on",
            "use_basis",
            "transformation",
        ):
            if not getattr(self, field).strip():
                raise ValueError(f"source {field} must be nonblank")
        expected = _source_fields(self.source_id)
        if (
            self.citation,
            self.locator,
            self.url,
            self.transformation,
            self.data_sha256,
            self.packaged_data_sha256,
        ) != expected:
            raise ValueError("source fields must match the exact retained source identity")
        if self.retrieved_on != SOURCE_RETRIEVED_ON or self.use_basis != SOURCE_USE_BASIS:
            raise ValueError("source provenance must match the exact retained source identity")
        if self.units != SOURCE_UNITS:
            raise ValueError("source units must match the exact retained table")


METHANE_SOURCE_V1 = SourceIdentity(
    source_id=METHANE_SOURCE_ID,
    citation=METHANE_SOURCE_CITATION,
    locator=METHANE_SOURCE_LOCATOR,
    url=METHANE_SOURCE_URL,
    retrieved_on=SOURCE_RETRIEVED_ON,
    use_basis=SOURCE_USE_BASIS,
    transformation=METHANE_SOURCE_TRANSFORMATION,
    units=SOURCE_UNITS,
    data_sha256=METHANE_DATA_SHA256,
    packaged_data_sha256=METHANE_PACKAGED_DATA_SHA256,
)
ETHANE_SOURCE_V1 = SourceIdentity(
    source_id=ETHANE_SOURCE_ID,
    citation=ETHANE_SOURCE_CITATION,
    locator=ETHANE_SOURCE_LOCATOR,
    url=ETHANE_SOURCE_URL,
    retrieved_on=SOURCE_RETRIEVED_ON,
    use_basis=SOURCE_USE_BASIS,
    transformation=ETHANE_SOURCE_TRANSFORMATION,
    units=SOURCE_UNITS,
    data_sha256=ETHANE_DATA_SHA256,
    packaged_data_sha256=ETHANE_PACKAGED_DATA_SHA256,
)


@dataclass(frozen=True, slots=True)
class SaturationObservation:
    row_id: str
    component_id: str
    temperature_k: float
    pressure_pa: float
    liquid_density_kg_m3: float
    source_id: str

    def __post_init__(self) -> None:
        if not self.row_id.strip():
            raise ValueError("row_id must be nonblank")
        if self.component_id not in ("methane", "ethane"):
            raise ValueError("component_id must be 'methane' or 'ethane'")
        _positive_finite(self.temperature_k, "temperature_k")
        _positive_finite(self.pressure_pa, "pressure_pa")
        _positive_finite(self.liquid_density_kg_m3, "liquid_density_kg_m3")
        expected_source = METHANE_SOURCE_ID if self.component_id == "methane" else ETHANE_SOURCE_ID
        if self.source_id != expected_source:
            raise ValueError("source_id does not match the pure component")


@dataclass(frozen=True, slots=True)
class PureSaturationDataset:
    dataset_id: str
    component_id: str
    temperature_unit: str
    pressure_unit: str
    liquid_density_unit: str
    source: SourceIdentity
    rows: tuple[SaturationObservation, ...]
    training_temperatures_k: tuple[float, ...]
    held_out_temperatures_k: tuple[float, ...]
    stress_temperatures_k: tuple[float, ...]

    def __post_init__(self) -> None:
        if self.component_id == "methane":
            expected = (
                "nist-webbook-methane-saturation-100-180-k-v1",
                METHANE_SOURCE_ID,
                METHANE_TEMPERATURES_K,
                METHANE_TRAINING_TEMPERATURES_K,
                METHANE_HELD_OUT_TEMPERATURES_K,
                (),
            )
        elif self.component_id == "ethane":
            expected = (
                "nist-webbook-ethane-saturation-100-280-k-v1",
                ETHANE_SOURCE_ID,
                ETHANE_TEMPERATURES_K,
                ETHANE_TRAINING_TEMPERATURES_K,
                ETHANE_HELD_OUT_TEMPERATURES_K,
                ETHANE_STRESS_TEMPERATURES_K,
            )
        else:
            raise ValueError("component_id must be 'methane' or 'ethane'")
        dataset_id, source_id, temperatures, training, held_out, stress = expected
        if self.dataset_id != dataset_id:
            raise ValueError("dataset_id does not match the admitted component table")
        if self.source.source_id != source_id:
            raise ValueError("dataset source does not match its component")
        if (self.temperature_unit, self.pressure_unit, self.liquid_density_unit) != (
            "K",
            "Pa",
            "kg/m3",
        ):
            raise ValueError("dataset units must be K, Pa, and kg/m3")
        row_ids = tuple(row.row_id for row in self.rows)
        if len(set(row_ids)) != len(row_ids):
            raise ValueError("dataset contains a duplicate row_id")
        observed_temperatures = tuple(row.temperature_k for row in self.rows)
        if any(
            right <= left
            for left, right in zip(observed_temperatures[:-1], observed_temperatures[1:], strict=True)
        ):
            raise ValueError("dataset temperatures must be strictly increasing")
        if observed_temperatures != temperatures:
            raise ValueError("dataset temperatures do not match the retained reporting grid")
        expected_row_ids = tuple(
            f"nist-{self.component_id}-sat-{int(temperature)}-k" for temperature in temperatures
        )
        if row_ids != expected_row_ids:
            raise ValueError("dataset row IDs must match their retained temperatures")
        if any(
            row.component_id != self.component_id or row.source_id != source_id
            for row in self.rows
        ):
            raise ValueError("every row must share the dataset component and source identity")
        if (
            self.training_temperatures_k,
            self.held_out_temperatures_k,
            self.stress_temperatures_k,
        ) != (training, held_out, stress):
            raise ValueError("dataset row partition does not match the admitted specification")
        if set(training) & set(held_out) or set(training) & set(stress) or set(held_out) & set(stress):
            raise ValueError("dataset row partitions must be disjoint")
        if set(training) | set(held_out) | set(stress) != set(temperatures):
            raise ValueError("dataset row partitions must cover every retained row")

    def _rows_at(self, temperatures: tuple[float, ...]) -> tuple[SaturationObservation, ...]:
        selected = frozenset(temperatures)
        return tuple(row for row in self.rows if row.temperature_k in selected)

    @property
    def training_rows(self) -> tuple[SaturationObservation, ...]:
        return self._rows_at(self.training_temperatures_k)

    @property
    def held_out_rows(self) -> tuple[SaturationObservation, ...]:
        return self._rows_at(self.held_out_temperatures_k)

    @property
    def stress_rows(self) -> tuple[SaturationObservation, ...]:
        return self._rows_at(self.stress_temperatures_k)

    @property
    def training_row_ids(self) -> tuple[str, ...]:
        return tuple(row.row_id for row in self.training_rows)


@dataclass(frozen=True, slots=True)
class PureSaturationFitSpecification:
    specification_id: str
    component_id: str
    dataset_id: str
    source_id: str
    expected_provider_fingerprint: str
    parameter_names: tuple[str, str, str]
    parameter_units: tuple[str, str, str]
    start: tuple[float, float, float]
    lower_bounds: tuple[float, float, float]
    upper_bounds: tuple[float, float, float]
    parameter_scales: tuple[float, float, float]
    fixed_amount_mol: float
    molar_mass_kg_per_mol: float
    residual_names: tuple[str, str, str, str]
    residual_weights: tuple[float, float, float, float]
    liquid_volume_bounds_m3: tuple[float, float]
    vapor_volume_bounds_m3: tuple[float, float]
    training_temperatures_k: tuple[float, float, float, float]
    max_num_iterations: int
    function_tolerance: float
    gradient_tolerance: float
    parameter_tolerance: float
    topology_relative_separation_min: float
    reporting_pressure_bounds_pa: tuple[float, float]
    confirmation_liquid_volume_start_multiplier: float
    confirmation_vapor_volume_start_multiplier: float
    confirmation_parameter_scaled_max_delta: float
    confirmation_cost_relative_delta: float
    reporting_pressure_scaled_residual_max: float
    reporting_chemical_potential_residual_max: float
    ceres_linear_solver: str
    ceres_num_threads: int
    ceres_logging: str

    def __post_init__(self) -> None:
        if self.component_id == "methane":
            expected_identity = (
                "pure-methane-saturation-lifted-volumes-v1",
                "nist-webbook-methane-saturation-100-180-k-v1",
                METHANE_SOURCE_ID,
                "sha256:5f836aa84935df70be2e5cffae51b178a7b797c2cee036e9ff47d8097ca94bbf",
                (1.08, 3.555744, 157.5315),
                0.016043,
                (1.5e-4, 0.1),
                METHANE_TRAINING_TEMPERATURES_K,
                (1.0e3, 1.0e7),
            )
        elif self.component_id == "ethane":
            expected_identity = (
                "pure-ethane-saturation-lifted-volumes-v1",
                "nist-webbook-ethane-saturation-100-280-k-v1",
                ETHANE_SOURCE_ID,
                "sha256:288fbcaa1304881c16f64c3a784eeed19b75c58cca4558f92a21268e5e91258a",
                (1.6069, 3.5206, 191.42),
                0.030070,
                (1.5e-4, 100.0),
                ETHANE_TRAINING_TEMPERATURES_K,
                (1.0, 1.0e7),
            )
        else:
            raise ValueError("component_id must be 'methane' or 'ethane'")
        if (
            self.specification_id,
            self.dataset_id,
            self.source_id,
            self.expected_provider_fingerprint,
            self.start,
            self.molar_mass_kg_per_mol,
            self.vapor_volume_bounds_m3,
            self.training_temperatures_k,
            self.reporting_pressure_bounds_pa,
        ) != expected_identity:
            raise ValueError("fit specification does not match the admitted component contract")
        if self.parameter_names != (
            "segment_count",
            "segment_diameter_angstrom",
            "dispersion_energy_over_k_kelvin",
        ) or self.parameter_units != ("1", "angstrom", "K"):
            raise ValueError("parameter names and units must match the provider coordinate contract")
        if any(
            not math.isfinite(value)
            for group in (self.start, self.lower_bounds, self.upper_bounds, self.parameter_scales)
            for value in group
        ):
            raise ValueError("parameter contract values must be finite")
        if any(
            not lower < start < upper
            for lower, start, upper in zip(
                self.lower_bounds, self.start, self.upper_bounds, strict=True
            )
        ):
            raise ValueError("every parameter start must lie strictly inside its bounds")
        if self.lower_bounds != (0.5, 2.0, 50.0) or self.upper_bounds != (3.5, 5.0, 400.0):
            raise ValueError("parameter bounds do not match the pure-saturation contract")
        if self.parameter_scales != (0.1, 0.1, 10.0):
            raise ValueError("parameter scales do not match the pure-saturation contract")
        _positive_finite(self.fixed_amount_mol, "fixed_amount_mol")
        if self.fixed_amount_mol != 1.0:
            raise ValueError("the pure-saturation slice fixes amount at exactly one mole")
        _positive_finite(self.molar_mass_kg_per_mol, "molar mass")
        if self.residual_names != (
            "liquid_pressure",
            "vapor_pressure",
            "chemical_potential_equality",
            "liquid_density",
        ):
            raise ValueError("residual ordering must match the lifted-volume formulation")
        if self.residual_weights != (0.25, 0.25, 0.25, 0.25):
            raise ValueError("residual weights must be the declared equal row normalization")
        if self.liquid_volume_bounds_m3 != (2.0e-5, 1.0e-4):
            raise ValueError("liquid volume bounds do not match the pure-saturation contract")
        if self.liquid_volume_bounds_m3[1] >= self.vapor_volume_bounds_m3[0]:
            raise ValueError("phase volume bounds must enforce liquid below vapor")
        if self.max_num_iterations != 500:
            raise ValueError("max_num_iterations does not match the pure-saturation contract")
        for value in (
            self.function_tolerance,
            self.gradient_tolerance,
            self.parameter_tolerance,
            self.topology_relative_separation_min,
            self.confirmation_liquid_volume_start_multiplier,
            self.confirmation_vapor_volume_start_multiplier,
            self.confirmation_parameter_scaled_max_delta,
            self.confirmation_cost_relative_delta,
            self.reporting_pressure_scaled_residual_max,
            self.reporting_chemical_potential_residual_max,
        ):
            _positive_finite(value, "solver tolerance")
        if (
            self.function_tolerance,
            self.gradient_tolerance,
            self.parameter_tolerance,
            self.topology_relative_separation_min,
        ) != (1.0e-10, 1.0e-10, 1.0e-10, 1.0e-3):
            raise ValueError("solver tolerances do not match the pure-saturation contract")
        if (
            self.confirmation_liquid_volume_start_multiplier,
            self.confirmation_vapor_volume_start_multiplier,
        ) != (1.01, 0.98):
            raise ValueError("confirmation start multipliers do not match the contract")
        if (
            self.confirmation_parameter_scaled_max_delta,
            self.confirmation_cost_relative_delta,
        ) != (1.0e-5, 1.0e-8):
            raise ValueError("confirmation acceptance thresholds do not match the contract")
        if (
            self.reporting_pressure_scaled_residual_max,
            self.reporting_chemical_potential_residual_max,
        ) != (1.0e-8, 1.0e-8):
            raise ValueError("reporting closure thresholds do not match the contract")
        if (self.ceres_linear_solver, self.ceres_num_threads, self.ceres_logging) != (
            "DENSE_QR",
            1,
            "SILENT",
        ):
            raise ValueError("Ceres execution controls do not match the contract")


def _fit_specification(
    *,
    component_id: str,
    specification_id: str,
    dataset_id: str,
    source_id: str,
    expected_provider_fingerprint: str,
    start: tuple[float, float, float],
    molar_mass_kg_per_mol: float,
    vapor_volume_bounds_m3: tuple[float, float],
    training_temperatures_k: tuple[float, float, float, float],
    reporting_pressure_bounds_pa: tuple[float, float],
) -> PureSaturationFitSpecification:
    return PureSaturationFitSpecification(
        specification_id=specification_id,
        component_id=component_id,
        dataset_id=dataset_id,
        source_id=source_id,
        expected_provider_fingerprint=expected_provider_fingerprint,
        parameter_names=(
            "segment_count",
            "segment_diameter_angstrom",
            "dispersion_energy_over_k_kelvin",
        ),
        parameter_units=("1", "angstrom", "K"),
        start=start,
        lower_bounds=(0.5, 2.0, 50.0),
        upper_bounds=(3.5, 5.0, 400.0),
        parameter_scales=(0.1, 0.1, 10.0),
        fixed_amount_mol=1.0,
        molar_mass_kg_per_mol=molar_mass_kg_per_mol,
        residual_names=(
            "liquid_pressure",
            "vapor_pressure",
            "chemical_potential_equality",
            "liquid_density",
        ),
        residual_weights=(0.25, 0.25, 0.25, 0.25),
        liquid_volume_bounds_m3=(2.0e-5, 1.0e-4),
        vapor_volume_bounds_m3=vapor_volume_bounds_m3,
        training_temperatures_k=training_temperatures_k,
        max_num_iterations=500,
        function_tolerance=1.0e-10,
        gradient_tolerance=1.0e-10,
        parameter_tolerance=1.0e-10,
        topology_relative_separation_min=1.0e-3,
        reporting_pressure_bounds_pa=reporting_pressure_bounds_pa,
        confirmation_liquid_volume_start_multiplier=1.01,
        confirmation_vapor_volume_start_multiplier=0.98,
        confirmation_parameter_scaled_max_delta=1.0e-5,
        confirmation_cost_relative_delta=1.0e-8,
        reporting_pressure_scaled_residual_max=1.0e-8,
        reporting_chemical_potential_residual_max=1.0e-8,
        ceres_linear_solver="DENSE_QR",
        ceres_num_threads=1,
        ceres_logging="SILENT",
    )


METHANE_SATURATION_FIT_V1 = _fit_specification(
    component_id="methane",
    specification_id="pure-methane-saturation-lifted-volumes-v1",
    dataset_id="nist-webbook-methane-saturation-100-180-k-v1",
    source_id=METHANE_SOURCE_ID,
    expected_provider_fingerprint=(
        "sha256:5f836aa84935df70be2e5cffae51b178a7b797c2cee036e9ff47d8097ca94bbf"
    ),
    start=(1.08, 3.555744, 157.5315),
    molar_mass_kg_per_mol=0.016043,
    vapor_volume_bounds_m3=(1.5e-4, 0.1),
    training_temperatures_k=METHANE_TRAINING_TEMPERATURES_K,
    reporting_pressure_bounds_pa=(1.0e3, 1.0e7),
)
ETHANE_SATURATION_FIT_V1 = _fit_specification(
    component_id="ethane",
    specification_id="pure-ethane-saturation-lifted-volumes-v1",
    dataset_id="nist-webbook-ethane-saturation-100-280-k-v1",
    source_id=ETHANE_SOURCE_ID,
    expected_provider_fingerprint=(
        "sha256:288fbcaa1304881c16f64c3a784eeed19b75c58cca4558f92a21268e5e91258a"
    ),
    start=(1.6069, 3.5206, 191.42),
    molar_mass_kg_per_mol=0.030070,
    vapor_volume_bounds_m3=(1.5e-4, 100.0),
    training_temperatures_k=ETHANE_TRAINING_TEMPERATURES_K,
    reporting_pressure_bounds_pa=(1.0, 1.0e7),
)


def _load_dataset(
    filename: str,
    component_id: str,
    display_name: str,
    source: SourceIdentity,
    dataset_id: str,
    temperatures: tuple[float, ...],
    training: tuple[float, ...],
    held_out: tuple[float, ...],
    stress: tuple[float, ...],
) -> PureSaturationDataset:
    data = files("epcsaft_regression").joinpath(f"data/{filename}").read_bytes()
    if hashlib.sha256(data).hexdigest() != source.packaged_data_sha256:
        raise ValueError(f"packaged {component_id} data SHA-256 does not match the source record")
    parsed = csv.reader(io.StringIO(data.decode("utf-8"), newline=""))
    if tuple(next(parsed)) != EXPECTED_HEADER:
        raise ValueError(f"{component_id} source data header or units changed")
    rows: list[SaturationObservation] = []
    for values in parsed:
        if len(values) != len(EXPECTED_HEADER):
            raise ValueError(f"{component_id} source row has a missing field")
        species, temperature, pressure, density, source_url = values
        if species != display_name or source_url != source.url:
            raise ValueError(f"{component_id} source identity changed")
        temperature_k = float(temperature)
        rows.append(
            SaturationObservation(
                row_id=f"nist-{component_id}-sat-{int(temperature_k)}-k",
                component_id=component_id,
                temperature_k=temperature_k,
                pressure_pa=float(pressure),
                liquid_density_kg_m3=float(density),
                source_id=source.source_id,
            )
        )
    if tuple(row.temperature_k for row in rows) != temperatures:
        raise ValueError(f"{component_id} source rows changed")
    return PureSaturationDataset(
        dataset_id=dataset_id,
        component_id=component_id,
        temperature_unit="K",
        pressure_unit="Pa",
        liquid_density_unit="kg/m3",
        source=source,
        rows=tuple(rows),
        training_temperatures_k=training,
        held_out_temperatures_k=held_out,
        stress_temperatures_k=stress,
    )


def load_pure_saturation_dataset(component_id: str) -> PureSaturationDataset:
    if type(component_id) is not str:
        raise TypeError("component_id must be an exact string")
    if component_id == "methane":
        return _load_dataset(
            "methane_saturation.csv",
            "methane",
            "Methane",
            METHANE_SOURCE_V1,
            "nist-webbook-methane-saturation-100-180-k-v1",
            METHANE_TEMPERATURES_K,
            METHANE_TRAINING_TEMPERATURES_K,
            METHANE_HELD_OUT_TEMPERATURES_K,
            (),
        )
    if component_id == "ethane":
        return _load_dataset(
            "ethane_saturation.csv",
            "ethane",
            "Ethane",
            ETHANE_SOURCE_V1,
            "nist-webbook-ethane-saturation-100-280-k-v1",
            ETHANE_TEMPERATURES_K,
            ETHANE_TRAINING_TEMPERATURES_K,
            ETHANE_HELD_OUT_TEMPERATURES_K,
            ETHANE_STRESS_TEMPERATURES_K,
        )
    raise ValueError("component_id must be 'methane' or 'ethane'")
