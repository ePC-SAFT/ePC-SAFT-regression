from __future__ import annotations

from dataclasses import dataclass
import csv
import hashlib
from importlib.resources import files
import io
import math


EXPECTED_DATA_SHA256 = "a5e16df3bf8ec78483fc340782cddc89ab8b284a9f6dfaecd6cda3ffde579227"
EXPECTED_PACKAGED_DATA_SHA256 = (
    "dec64d5a6cac414a4a92393a0d728fa27c02135c6a159d0d1881d7b6dde6d26c"
)
SOURCE_ID = "nist-webbook-srd69-methane-saturation"
SOURCE_URL = (
    "https://webbook.nist.gov/cgi/fluid.cgi?Action=Data&Wide=on&ID=C74828&"
    "Type=SatP&Digits=8&THigh=180&TLow=100&TInc=10&RefState=DEF&TUnit=K&"
    "PUnit=Pa&DUnit=kg%2Fm3&HUnit=kJ%2Fmol&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm"
)
SOURCE_CITATION = "NIST Chemistry WebBook, SRD 69, methane (CAS 74-82-8) fluid properties"
SOURCE_LOCATOR = "Saturation properties query, 100 K through 180 K in 10 K increments"
SOURCE_RETRIEVED_ON = "2026-07-17"
SOURCE_USE_BASIS = (
    "NIST Standard Reference Data retained as compact source-backed candidate evidence; "
    "redistribution and use remain subject to the NIST SRD terms"
)
SOURCE_TRANSFORMATION = (
    "Exact retained CSV fields and decimal strings; CRLF line endings normalized to LF. "
    f"Retained source SHA-256: {EXPECTED_DATA_SHA256}; packaged SHA-256: "
    f"{EXPECTED_PACKAGED_DATA_SHA256}."
)
SOURCE_UNITS = (
    ("temperature", "K"),
    ("pressure", "Pa"),
    ("saturated_liquid_mass_density", "kg/m3"),
)
EXPECTED_HEADER = ("species", "T_K", "p_sat_Pa", "rho_sat_liq_kg_m3", "source")
EXPECTED_TEMPERATURES = tuple(float(value) for value in range(100, 181, 10))
TRAINING_TEMPERATURES = (110.0, 130.0, 150.0, 170.0)


def _positive_finite(value: float, field: str) -> None:
    if not math.isfinite(value):
        raise ValueError(f"{field} must be finite")
    if value <= 0.0:
        raise ValueError(f"{field} must be positive")


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
        if (
            self.source_id,
            self.citation,
            self.locator,
            self.url,
            self.retrieved_on,
            self.use_basis,
            self.transformation,
        ) != (
            SOURCE_ID,
            SOURCE_CITATION,
            SOURCE_LOCATOR,
            SOURCE_URL,
            SOURCE_RETRIEVED_ON,
            SOURCE_USE_BASIS,
            SOURCE_TRANSFORMATION,
        ):
            raise ValueError("source fields must match the exact retained source identity")
        if self.units != SOURCE_UNITS:
            raise ValueError("source units must match the exact retained table")
        if self.data_sha256 != EXPECTED_DATA_SHA256:
            raise ValueError("source data SHA-256 does not match the admitted retained table")


@dataclass(frozen=True, slots=True)
class SaturationObservation:
    row_id: str
    species: str
    temperature_k: float
    pressure_pa: float
    liquid_density_kg_m3: float
    source_id: str

    def __post_init__(self) -> None:
        if not self.row_id.strip():
            raise ValueError("row_id must be nonblank")
        if self.species != "methane":
            raise ValueError("the first slice supports methane only")
        _positive_finite(self.temperature_k, "temperature_k")
        _positive_finite(self.pressure_pa, "pressure_pa")
        _positive_finite(self.liquid_density_kg_m3, "liquid_density_kg_m3")
        if self.source_id != SOURCE_ID:
            raise ValueError("source_id must identify the retained NIST methane table")


@dataclass(frozen=True, slots=True)
class MethaneSaturationDataset:
    dataset_id: str
    species: str
    temperature_unit: str
    pressure_unit: str
    liquid_density_unit: str
    source: SourceIdentity
    rows: tuple[SaturationObservation, ...]
    training_row_ids: tuple[str, ...]

    def __post_init__(self) -> None:
        if self.dataset_id != "nist-webbook-methane-saturation-100-180-k-v1":
            raise ValueError("dataset_id must identify the first admitted methane table")
        if self.species != "methane":
            raise ValueError("dataset species must be methane")
        if (self.temperature_unit, self.pressure_unit, self.liquid_density_unit) != (
            "K",
            "Pa",
            "kg/m3",
        ):
            raise ValueError("dataset units must be K, Pa, and kg/m3")
        row_ids = tuple(row.row_id for row in self.rows)
        if len(set(row_ids)) != len(row_ids):
            raise ValueError("dataset contains a duplicate row_id")
        if len(self.rows) != len(EXPECTED_TEMPERATURES):
            raise ValueError("dataset must contain exactly nine reporting rows")
        temperatures = tuple(row.temperature_k for row in self.rows)
        if any(right <= left for left, right in zip(temperatures[:-1], temperatures[1:], strict=True)):
            raise ValueError("dataset temperatures must be strictly increasing")
        if temperatures != EXPECTED_TEMPERATURES:
            raise ValueError("dataset temperatures do not match the retained reporting grid")
        expected_row_ids = tuple(
            f"nist-methane-sat-{int(temperature)}-k" for temperature in EXPECTED_TEMPERATURES
        )
        if row_ids != expected_row_ids:
            raise ValueError("dataset row IDs must match their retained temperatures")
        if any(row.source_id != self.source.source_id for row in self.rows):
            raise ValueError("every row must share the dataset source identity")
        expected_training = tuple(
            f"nist-methane-sat-{int(temperature)}-k" for temperature in TRAINING_TEMPERATURES
        )
        if self.training_row_ids != expected_training:
            raise ValueError("dataset training partition must be 110, 130, 150, and 170 K")
        if any(row_id not in row_ids for row_id in self.training_row_ids):
            raise ValueError("dataset training partition references an unknown row")

    @property
    def training_rows(self) -> tuple[SaturationObservation, ...]:
        selected = frozenset(self.training_row_ids)
        return tuple(row for row in self.rows if row.row_id in selected)


@dataclass(frozen=True, slots=True)
class MethaneFitSpecification:
    specification_id: str
    parameter_names: tuple[str, str, str]
    parameter_units: tuple[str, str, str]
    start: tuple[float, float, float]
    lower_bounds: tuple[float, float, float]
    upper_bounds: tuple[float, float, float]
    parameter_scales: tuple[float, float, float]
    fixed_amount_mol: float
    methane_molar_mass_kg_per_mol: float
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
        if self.specification_id != "pure-methane-saturation-lifted-volumes-v1":
            raise ValueError("unsupported methane fit specification")
        if self.parameter_names != (
            "segment_count",
            "segment_diameter_angstrom",
            "dispersion_energy_over_k_kelvin",
        ) or self.parameter_units != ("1", "angstrom", "K"):
            raise ValueError("parameter names and units must match the provider coordinate contract")
        if any(not math.isfinite(value) for group in (
            self.start,
            self.lower_bounds,
            self.upper_bounds,
            self.parameter_scales,
        ) for value in group):
            raise ValueError("parameter contract values must be finite")
        if any(not lower < start < upper for lower, start, upper in zip(
            self.lower_bounds, self.start, self.upper_bounds, strict=True
        )):
            raise ValueError("every parameter start must lie strictly inside its bounds")
        if any(scale <= 0.0 for scale in self.parameter_scales):
            raise ValueError("parameter scales must be positive")
        _positive_finite(self.fixed_amount_mol, "fixed_amount_mol")
        if self.fixed_amount_mol != 1.0:
            raise ValueError("the first slice fixes amount at exactly one mole")
        _positive_finite(self.methane_molar_mass_kg_per_mol, "methane molar mass")
        if self.residual_names != (
            "liquid_pressure",
            "vapor_pressure",
            "chemical_potential_equality",
            "liquid_density",
        ):
            raise ValueError("residual ordering must match the first lifted-volume formulation")
        if self.residual_weights != (0.25, 0.25, 0.25, 0.25):
            raise ValueError("residual weights must be the declared equal row normalization")
        if self.liquid_volume_bounds_m3 != (2.0e-5, 1.0e-4):
            raise ValueError("liquid volume bounds do not match the first slice")
        if self.vapor_volume_bounds_m3 != (1.5e-4, 0.1):
            raise ValueError("vapor volume bounds do not match the first slice")
        if self.liquid_volume_bounds_m3[1] >= self.vapor_volume_bounds_m3[0]:
            raise ValueError("phase volume bounds must enforce liquid below vapor")
        if self.training_temperatures_k != TRAINING_TEMPERATURES:
            raise ValueError("fit specification training temperatures are unsupported")
        if self.max_num_iterations <= 0:
            raise ValueError("max_num_iterations must be positive")
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
        if self.reporting_pressure_bounds_pa != (1.0e3, 1.0e7):
            raise ValueError("reporting pressure bounds do not match the first slice")
        if (
            self.confirmation_liquid_volume_start_multiplier,
            self.confirmation_vapor_volume_start_multiplier,
        ) != (1.01, 0.98):
            raise ValueError("confirmation start multipliers do not match the first slice")
        if (
            self.confirmation_parameter_scaled_max_delta,
            self.confirmation_cost_relative_delta,
        ) != (1.0e-5, 1.0e-8):
            raise ValueError("confirmation acceptance thresholds do not match the first slice")
        if (
            self.reporting_pressure_scaled_residual_max,
            self.reporting_chemical_potential_residual_max,
        ) != (1.0e-8, 1.0e-8):
            raise ValueError("reporting closure thresholds do not match the first slice")
        if (
            self.ceres_linear_solver,
            self.ceres_num_threads,
            self.ceres_logging,
        ) != ("DENSE_QR", 1, "SILENT"):
            raise ValueError("Ceres execution controls do not match the first slice")


METHANE_FIT_SPECIFICATION_V1 = MethaneFitSpecification(
    specification_id="pure-methane-saturation-lifted-volumes-v1",
    parameter_names=(
        "segment_count",
        "segment_diameter_angstrom",
        "dispersion_energy_over_k_kelvin",
    ),
    parameter_units=("1", "angstrom", "K"),
    start=(1.08, 3.555744, 157.5315),
    lower_bounds=(0.5, 2.0, 50.0),
    upper_bounds=(3.5, 5.0, 400.0),
    parameter_scales=(0.1, 0.1, 10.0),
    fixed_amount_mol=1.0,
    methane_molar_mass_kg_per_mol=0.016043,
    residual_names=(
        "liquid_pressure",
        "vapor_pressure",
        "chemical_potential_equality",
        "liquid_density",
    ),
    residual_weights=(0.25, 0.25, 0.25, 0.25),
    liquid_volume_bounds_m3=(2.0e-5, 1.0e-4),
    vapor_volume_bounds_m3=(1.5e-4, 0.1),
    training_temperatures_k=TRAINING_TEMPERATURES,
    max_num_iterations=500,
    function_tolerance=1.0e-10,
    gradient_tolerance=1.0e-10,
    parameter_tolerance=1.0e-10,
    topology_relative_separation_min=1.0e-3,
    reporting_pressure_bounds_pa=(1.0e3, 1.0e7),
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


def load_methane_dataset() -> MethaneSaturationDataset:
    data = files("epcsaft_regression").joinpath("data/methane_saturation.csv").read_bytes()
    if hashlib.sha256(data).hexdigest() != EXPECTED_PACKAGED_DATA_SHA256:
        raise ValueError("packaged methane data SHA-256 does not match the distilled evidence")
    parsed = csv.reader(io.StringIO(data.decode("utf-8"), newline=""))
    header = tuple(next(parsed))
    if header != EXPECTED_HEADER:
        raise ValueError("methane source data header or units changed")
    rows: list[SaturationObservation] = []
    for values in parsed:
        if len(values) != len(EXPECTED_HEADER):
            raise ValueError("methane source row has a missing field")
        species, temperature, pressure, density, source_url = values
        if species != "Methane" or source_url != SOURCE_URL:
            raise ValueError("methane source identity changed")
        temperature_k = float(temperature)
        rows.append(
            SaturationObservation(
                row_id=f"nist-methane-sat-{int(temperature_k)}-k",
                species="methane",
                temperature_k=temperature_k,
                pressure_pa=float(pressure),
                liquid_density_kg_m3=float(density),
                source_id=SOURCE_ID,
            )
        )
    source = SourceIdentity(
        source_id=SOURCE_ID,
        citation=SOURCE_CITATION,
        locator=SOURCE_LOCATOR,
        url=SOURCE_URL,
        retrieved_on=SOURCE_RETRIEVED_ON,
        use_basis=SOURCE_USE_BASIS,
        transformation=SOURCE_TRANSFORMATION,
        units=SOURCE_UNITS,
        data_sha256=EXPECTED_DATA_SHA256,
    )
    training_row_ids = tuple(
        f"nist-methane-sat-{int(temperature)}-k" for temperature in TRAINING_TEMPERATURES
    )
    return MethaneSaturationDataset(
        dataset_id="nist-webbook-methane-saturation-100-180-k-v1",
        species="methane",
        temperature_unit="K",
        pressure_unit="Pa",
        liquid_density_unit="kg/m3",
        source=source,
        rows=tuple(rows),
        training_row_ids=training_row_ids,
    )
