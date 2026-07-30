# General Regression quickstart

The public direct-EOS path is `prepare_fit -> preflight -> fit -> export`.
Preparation resolves capability IDs, EOS units, parameter/topology
fingerprints, and transformed-row hashes from the installed model. The caller
still supplies every scientific choice: source, rows, objective, bounds,
affine scales, complete starts, partitions, and lifted-volume contracts.

```python
from epcsaft import Mixture, Parameters
from epcsaft_regression import (
    AcquisitionClass, AffineParameterTransform, ConfirmationControls,
    FixedCompositionVleObservation, ObjectiveContract, ObservationDataset,
    PairParameterIdentity, ParameterFamily, ParameterRequest, RankControls,
    ResultContext, RowProvenance, SolverControls, SourceInput, prepare_fit,
    support_view,
)

model = Mixture(Parameters.from_catalog(
    "gross-2001-methane-ethane",
    components=("methane", "ethane"), version=1,
))
assert any(
    capability.installed_ready
    and getattr(capability, "family", None) is ParameterFamily.K_IJ
    for capability in support_view(model)
)
records = ({
    "row_id": "may2015-ch4-c2h6-002",
    "source_id": "may-2015",
    "source_locator": "May et al. (2015), retained row 002",
    "component_ids": ["methane", "ethane"],
    "temperature_k": 203.22,
    "pressure_pa": 2_124_000.0,
    "liquid_mole_fraction_first": 0.3653,
    "vapor_mole_fraction_first": 0.8667,
    "pressure_scale_pa": 2_124_000.0,
    "chemical_potential_scales": [1.0, 1.0],
    "liquid_volume_origin_m3_per_mol": 6.0e-5,
    "liquid_volume_start_m3_per_mol": 6.5e-5,
    "liquid_volume_bounds_m3_per_mol": [2.0e-5, 1.0e-4],
    "vapor_volume_origin_m3_per_mol": 9.0e-4,
    "vapor_volume_start_m3_per_mol": 1.0e-3,
    "vapor_volume_bounds_m3_per_mol": [1.0e-4, 1.0e-2],
    "partition": "training",
},)
dataset = ObservationDataset.from_records(
    FixedCompositionVleObservation, records,
    source=SourceInput(
        "may-2015", "May et al. (2015), methane/ethane VLE",
        "retained source artifact",
        "5cd1e74925a3c6504f5106dcf911f2cae2d6e99a5133fccc20454d8991bdbc7f",
        "No transformation for this retained row.",
        "T/K, P/Pa, mole fractions",
        "Regression derivative-contract evidence",
        "pressure by observed P; mu/RT dimensionless",
    ),
    objective=ObjectiveContract(
        "fixed_composition_vle", "native_scaled_least_squares",
        "scales embedded in each validated row",
        "independent rows; no covariance supplied", "squared", (), "fail",
    ),
    row_provenance={"may2015-ch4-c2h6-002": RowProvenance(
        AcquisitionClass.DIRECT_MEASUREMENT, "unique source row", "included",
        "outside critical-region exclusion", "not censored",
        "retained; no outlier rule applied",
    )},
)
prepared = prepare_fit(
    model, datasets=(dataset,),
    parameters=(ParameterRequest(
        ParameterFamily.K_IJ, PairParameterIdentity("methane", "ethane"),
        AffineParameterTransform(0.0, 0.01), -0.15, 0.10,
    ),),
    parameter_slot_indices=(0,),
    start_vectors=((0.0,), (-0.05,), (0.05,)),
    solver=SolverControls(50, 30.0, 1e-12, 1e-12, 1e-12),
    rank=RankControls(1e10),
    confirmation=ConfirmationControls(1e-5, 1e-8),
)
report = prepared.preflight()  # exact evaluations; no Ceres solve
if not report.ready:
    raise RuntimeError(report.reasons)
result = prepared.fit()
canonical_json = result.to_json_bytes(prepared=prepared, context=ResultContext())
```

Preflight reports `N`, `Q`, `R`, partitions, derivative completeness, full and
nuisance-projected local rank/conditioning, start-bound status, and failure
reasons. It does not report practical/global identifiability, uncertainty,
prediction, or acceptance. Repeated or insensitive row designs therefore fail
rank preflight before Ceres.

`support_view(model)` returns the installed EOS capability records directly,
including units, identity and observation contracts, derivative order,
`installed_ready`, and `unsupported_reason`. The human-readable family matrix
is in
[`general-parameter-regression.md`](science/general-parameter-regression.md#definitive-parameter-family-and-data-requirement-matrix).

## Parameter blocks and association

- Scalar fits use one request and one slot.
- Joint pure `(m, sigma, epsilon/k)` uses those families and slots in that
  order, complete three-value starts, and pure-saturation rows.
- Fixed neutral-pure-2B uses `(m, sigma, epsilon/k, epsilon_AB/k, kappa_AB)`
  in that order with pure vapor-pressure, fixed-pressure-density, and/or
  combined saturation rows. Topology is fixed EOS input. Issue #28 owns
  nuisance-reoptimized profile and accepted-region evidence. A constrained
  subset such as externally fixed `kappa_AB` must declare its source and
  validity domain; it is not universal.
- Direct-observable families use their MIAC, solvation-Gibbs, or relative-
  permittivity-ratio rows and no lifted variables.

Induced and resolved cross association are not current pure-fit coordinates.
Issue #34 owns future fixed alternatives after EOS advertises their exact
identities, transforms, topology fingerprints, and derivatives.

## Record templates

Mappings require every dataclass field and reject unknown fields.

| Observation | Scientific fields beyond IDs/source/partition |
|---|---|
| Fixed-composition VLE | ordered pair, `T`, `P`, liquid/vapor composition, pressure and chemical-potential scales, phase volume origins/starts/bounds |
| Pure saturation | component, `T`, `P_sat`, liquid density, molar mass, residual scales, phase volume origins/starts/bounds |
| Pure vapor pressure | component, `T`, `P_vap`, residual scales, separated phase volume origins/starts/bounds |
| Fixed-pressure pure density | component, `T`, `P`, density, molar mass, residual scales, volume origin/start/bounds |
| MIAC / aqueous pair MIAC | ordered model identity, active component or pair, fixed pair context where applicable, `T`, `P`, molality, target and scale |
| Solvation Gibbs / ion-pair solvation | ordered model identity, observed ion/active coordinate, fixed pair context where applicable, `T`, `P`, target and scale |
| Relative-permittivity ratio | solvent and ordered ternary identity, `T`, `P`, total-ion mole fraction, target and scale |

Every row requires acquisition class and explicit duplicate, exclusion,
critical-region, censoring, and outlier decisions. Correlation rows additionally
require equation, coefficients, units, validity interval, exact sampling grid,
and transformation. Dense correlation samples are not independent experiments.

Export returns data only. It writes no receipt, mutates no EOS catalog, and
grants no production, predictive, or scientific authority.
