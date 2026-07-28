# General ePC-SAFT Parameter Regression

Status: neutral-binary `k_ij`/`l_ij` and scalar pure
`m`/`sigma`/`epsilon-k` families fit-ready; broader families pending

Date: 2026-07-27

## Objective

Regression shall support source-bound fitting of any admitted continuous
ePC-SAFT parameter for any Provider-supported component set when the caller
supplies an informative dataset and a complete numerical contract.

Published methane, ethane, Figiel, and later literature workflows are
reference campaigns. They provide positive or negative evidence for an exact
parameter-family domain but do not restrict the runtime to those components,
mixtures, or papers. One campaign does not prove scientific accuracy for every
mixture.

This design supersedes the earlier assumption that every paper-specific
campaign must remain a separate regression surface. It does not claim that
every Provider parameter is fit-ready today.

## Ownership boundary

Provider remains the only owner of ePC-SAFT equations, resolved parameter
records, thermodynamic states, density closure, and nonlinear derivatives.
Regression owns:

- source-bound observation contracts;
- active parameter selection, bounds, scales, transforms, starts, and sharing;
- residual construction and Ceres execution;
- exact residual-Jacobian assembly from Provider derivatives;
- rank, conditioning, active-bound, confirmation, and row-accounting evidence;
- immutable fit results and fitted-versus-reference comparisons.

Regression must not accept arbitrary executable residual plugins, copied EOS
equations, a numerical production derivative backend, or a mutable parameter
registry. The supported vocabulary for any fit is a closed, versioned typed
schema. Provider may advertise a versioned superset. Regression accepts only
descriptor members its installed version understands, reports the remainder
as unsupported, and rejects any request that references an unknown member
rather than registering it dynamically.

## Capability meanings

The package shall distinguish these states for every parameter family:

1. `REPRESENTED`: Provider can resolve and evaluate a model containing the
   parameter.
2. `DERIVATIVE_READY`: an installed Provider artifact exposes the exact active
   derivative required by at least one observation contract.
3. `FIT_READY`: Regression can construct and solve the corresponding
   source-bound problem for caller-supplied components and rows.
4. `REFERENCE_VALIDATED`: at least one source-backed installed-artifact
   campaign has passed its declared gates and retained producer, reviewer,
   artifact, command, and evidence identities.
5. `ADMITTED`: an accepted promotion receipt binds the exact capability,
   producer and distinct reviewer, source, installed artifacts, evidence
   subject, and user approval. Documentation, a merge, or a converged fit
   cannot set this state.

An example dataset may advance validation evidence without narrowing the
component domain of a `FIT_READY` family.

## Typed parameter identity

An active coordinate is identified by one closed identity variant, not by a
paper or array position:

```text
ParameterCoordinate
  family: ParameterFamily
  identity:
    ComponentParameter(component_id)
    | PairParameter(component_id_a, component_id_b)
    | AssociationParameter((component_id_a, site_id_a),
                           (component_id_b, site_id_b))
    | ModelParameter(model_parameter_id)
    | CorrelationCoefficient(component_id, correlation_family,
                             coefficient_kind, term_index)
  capability_id
  provider_parameter_fingerprint
  provider_topology_fingerprint
  transform: Affine(origin, scale)       # v1; p = origin + scale*z
  starts
  lower_bound
  upper_bound
```

Each family permits exactly one identity variant. Pair endpoints and
association endpoints use one documented canonical lexical order; duplicate
canonical identities are rejected. Provider's capability descriptor supplies
the canonical physical unit and an explicit mapping from these identities to
its active coordinate order. A model, formulation, topology, transform, or
parameter-fingerprint mismatch fails before evaluation.

Bounds, affine scale, and primary/confirmation starts are mandatory source or
workflow inputs; the runtime supplies no chemistry defaults. Multiple
observations may share one coordinate, and multiple coordinates may be fitted
together. A future non-affine transform requires a new closed transform kind
and its exact `dp/dz`; v1 must not infer one from parameter sign or bounds.

Discrete model choices, charge number, component identity, site topology, and
electrolyte formulation are specifications, not continuous fit coordinates.

### Complete Provider-record inventory

The current Provider record vocabulary is the starting inventory, not proof of
derivative readiness:

| Provider record family | Regression treatment |
|---|---|
| `molar_mass` | fixed measured input; not an adjustable ePC-SAFT coordinate |
| `charge_number` | fixed discrete identity; not a continuous coordinate |
| `electrolyte_formulation` | discrete model selection; not fitted |
| `segment_count`, `segment_diameter`, `dispersion_energy_over_k` | continuous component coordinates |
| `relative_permittivity`, `born_diameter`, `solvation_factor` | continuous component coordinates when an exact observation derivative is advertised |
| `schreckenberg_dielectric_volume`, `schreckenberg_dielectric_temperature` | continuous component/correlation coordinates when exact derivatives are advertised |
| `zuber_ion_suppression_coefficient` | continuous component coordinate when exact derivatives are advertised |
| `rueben_dipole_scaling`, `rueben_polarizability_scaling`, `rueben_correlation_integral_parameter` | continuous polar-model coordinates when exact derivatives are advertised |
| `k_ij`, `l_ij` | continuous unordered component-pair coordinates |
| `association_energy_over_k`, `association_volume` | continuous association-endpoint coordinates |
| `dielectric_ion_suppression_coefficient`, `ionic_region_relative_permittivity` | continuous model coordinates when exact derivatives are advertised |
| correlation terms for `segment_diameter`, `relative_permittivity`, and `solvation_factor` | individually named coefficient coordinates; no whole-correlation opaque fit |

New Provider record families require an explicit schema revision in both
packages. Regression must never infer fit support merely because a record can
be resolved.

## Parameter and data map

The data column names below describe observation families, not guaranteed
identifiability. Every submitted problem must still pass numerical sensitivity
and rank gates.

| Parameter family | Identity | Candidate informative observations | Current state |
|---|---|---|---|
| `segment_count` (`m`) | component | pure saturation pressure/density, PVT density, phase-equilibrium and caloric observations | fit-ready one-family-at-a-time for caller-supplied pure-saturation pressure/liquid-density rows on an advertised neutral pure capability; other observations pending |
| `segment_diameter` (`sigma`) | component or declared correlation coefficient | liquid density/PVT, saturation density, phase equilibrium | fit-ready one-family-at-a-time for caller-supplied pure-saturation pressure/liquid-density rows on an advertised neutral pure capability; correlation coefficients and other observations pending |
| `dispersion_energy_over_k` (`epsilon/k`) | component | vapor pressure, phase equilibrium, caloric and PVT observations | fit-ready one-family-at-a-time for caller-supplied pure-saturation pressure/liquid-density rows on an advertised neutral pure capability; other observations pending |
| `k_ij` | unordered component pair | fixed-composition VLE/LLE, activity/fugacity coefficients, MIAC when the electrolyte formulation makes that pair active | fit-ready for caller-supplied fixed-composition VLE rows on any installed Provider model advertising the neutral, nonassociating binary capability; other observation/model domains remain pending |
| `l_ij` | unordered component pair | density, excess volume, and phase-equilibrium data sensitive to cross-size mixing | fit-ready for caller-supplied fixed-composition VLE rows on any installed Provider model advertising the neutral, nonassociating binary capability; other observation/model domains remain pending |
| `association_energy_over_k` | association endpoint pair | hydrogen-bond-sensitive VLE/LLE, solvation, enthalpy, and related phase-property observations | represented by Provider; exact active derivative and Regression surface pending |
| `association_volume` | association endpoint pair | hydrogen-bond-sensitive VLE/LLE, density, solvation, and caloric observations | represented by Provider; exact active derivative and Regression surface pending |
| `k_hb_ij` | source-defined cross-association combining-rule coordinate | the same cross-association-sensitive observations, with enough composition/temperature variation to separate it from pure association parameters | not a current Provider record family; requires an explicit active combining-rule coordinate or a declared transform to resolved cross-association energy |
| `born_diameter` | ion component | source-defined single-ion solvation Gibbs energy or another explicitly admitted ionic observable | five-ion fit demonstrated under one fixed `(T,P,reference-path)` convention; general request pending |
| `solvation_factor` | applicable component | MIAC, osmotic/activity, or solvation observations under a declared reference sequence | one NaBr water-factor fit demonstrated under one fixed reference path; general request pending |
| dielectric and ion-suppression coefficients | component/model/correlation term | dielectric, MIAC, osmotic/activity, solvation, or phase observations over a rank-sufficient state range | represented families vary; exact active derivatives and Regression surface pending |
| polar coefficients | component/correlation term | polar-mixture VLE, PVT, caloric, dielectric, or other source-defined polar observables | represented families vary; exact active derivatives and Regression surface pending |
| temperature-dependent coefficients | named correlation term | observations spanning enough temperatures to identify every active coefficient | represented correlations vary; exact active derivatives and Regression surface pending |

The table describes physically relevant evidence, not an automatic weighting
policy. Experimental uncertainty may define a weight only when the source
actually reports it and the request declares that interpretation.

### `k_hb_ij` rule

`k_hb_ij` must not be inferred from a resolved cross-association energy and
silently fitted under another name. A fit may use it only when Provider exposes
the source-defined combining rule, its component/site identity, a versioned
transform identifier and fingerprint, and its exact total derivative including
the transform chain rule. Otherwise the caller may fit the resolved
`association_energy_over_k` coordinate, and the result must retain that
different meaning.

## Source-bound observations

Dataset-level provenance is stored once:

```text
SourceDescriptor
  source_id
  citation and durable locator
  source artifact SHA-256
  canonical transformed-dataset SHA-256
  transformation record
  units and composition/property bases
  license or use basis
  tolerance and residual-scale rationale
```

Every observation row references that descriptor and contains:

```text
row_id
observation_type
component_ids and exact order
T and P specifications
composition and basis
phase or reference-state convention
observed value and unit
residual scale and its rationale
partition: TRAINING | HELD_OUT | STRESS
source_id and row locator
```

The canonical transformed-dataset hash binds every row to exactly one
partition. Row IDs are globally unique, partitions are disjoint, and only
`TRAINING` rows may create Ceres residual blocks. `HELD_OUT` and `STRESS` rows
are evaluated after the primary fit without refitting. Any overlap,
relabeling, or attempt to pass an untouched row into Ceres fails closed.
Without a predeclared acceptance cutoff, their errors remain descriptive and
predictive status remains not adjudicated.

Initial typed observation contracts are:

- pure saturation coexistence;
- single-phase scalar property or density;
- fixed-composition phase equilibrium with observed phase compositions and
  lifted phase state variables;
- activity, fugacity, osmotic, or mean-ionic-activity coefficient;
- solvation Gibbs energy.

Positive equality, linear aggregate, and one-sided censored observations may
be added as typed contracts when a real source requires them. The current
dependency doctrine permits Regression to consume Provider only. Therefore
predicted bubble/dew/flash, reactive, or other eliminated-equilibrium
observations are not authorized by this design. They require a future doctrine
amendment and separately accepted value/sensitivity transport contract;
Regression must not import Equilibrium, run a second equilibrium
implementation, or finite-difference a black-box solve.

The loader rejects missing units, ambiguous composition bases, duplicate row
identities, unsupported phases, nonfinite values, missing source identities,
nonfinite or nonpositive residual scales, scales dimensionally incompatible
with their observed values, source/hash mismatches, and
parameter/observation combinations that the installed artifact does not
advertise. In canonical value units, every scale must satisfy
`0 < s_q < infinity`.

## Exact derivative contract

Provider `1e571ab0a84603a51ed6994b14286f683fb12b88` supplies the first two
general model-bound capability descriptors for neutral binary `k_ij` and
`l_ij`. Provider `86983ff` adds three component-identity descriptors for a
neutral, nonassociating, constant-diameter pure model, backed by the exact
`(n,V,active_parameter)` value/gradient/Hessian callback. Regression selects
the requested known
descriptor from that finite Provider capability set, reports unknown
descriptors as unsupported, and rejects an unknown capability request. Later
families must supply the same closed metadata:

- schema and capability identifiers;
- parameter family and identity shape;
- accepted observation or phase-block contract;
- full state and active coordinate order and units;
- row-major tensor layout and exact buffer dimensions;
- identity-to-coordinate mapping;
- returned value units;
- exact derivative order;
- fixed topology and domain limits;
- local/global scope and explicit non-global claims;
- phase/branch stability and discovery status where applicable;
- certificate and tolerance-set identity;
- explicit unsupported-model, domain, topology, branch, derivative, and
  evaluation failure boundaries;
- Provider source, parameter, topology, transform, and artifact fingerprints;
- accepted-artifact status and receipt identity, or an explicit
  authority-neutral status.

Missing or unknown descriptor metadata prevents `DERIVATIVE_READY` and
`FIT_READY`; Regression does not fill it from local defaults.

The first descriptor is derivative-ready and authority-neutral. Regression's
typed fixed-composition VLE contract makes that exact family `FIT_READY`; it
does not make the capability reference-validated or admitted.

Direct observable residuals consume Provider values and first total
derivatives. For solver coordinates `z_k`, physical parameters `p_j`, and
residual scale `s_q`:

```text
r_q = (y_q(p) - y_q,observed) / s_q
J_qk = sum_j[(dy_q/dp_j) * (dp_j/dz_k)] / s_q
```

`dp/dz` owns affine scales and any declared parameter sharing. It is never
silently assumed diagonal.

Lifted phase-equilibrium residuals consume Provider phase-potential values,
gradients, and Hessians. The Hessian supplies exact derivatives of pressure
and chemical-potential residuals with respect to lifted state variables and
active parameters; the capability's coordinate order and row-major layout are
authoritative. No third derivative is required for a Ceres Jacobian whose
residuals are values or first gradients of that phase potential.

For a scalar pure-saturation row with lifted
`V_L = V_L,origin exp(u_L)` and `V_V = V_V,origin exp(u_V)`, the four
residuals are

```text
r_P,L   = (P_L - P_observed) / s_P
r_P,V   = (P_V - P_observed) / s_P
r_mu    = (Phi_n,L - Phi_n,V) / s_mu
r_rho,L = (M/V_L - rho_L,observed) / s_rho
```

The exact parameter column uses
`dP/dtheta = -RT Phi_(V,theta)` and
`d(Phi_n)/dtheta = Phi_(n,theta)`. The lifted-volume columns use the
corresponding `Phi_(V,V)` and `Phi_(n,V)` entries times `V`; the density
column is `-(M/V_L)/s_rho`. Thus the Provider Hessian is sufficient and no
third derivative or density-root sensitivity is introduced. Regression rejects
out-of-bound or inverted volumes and nonpositive mechanical-stability
curvature before accepting an evaluation.

Provider may use diagnosed density closure for a direct density/property
contract. Branch identity, closure residual, and conditioning then remain
Provider diagnostics returned with the observation. Regression does not
reimplement the closure.

Unsupported family/observation pairs fail before Ceres starts.

## One engine and result

The intended public surface is:

```text
fit_parameters(problem: RegressionProblem) -> RegressionResult
parameter_capabilities(model) -> tuple[ParameterCapability, ...]
```

`RegressionProblem` owns one ordered parameter block, typed observations,
source identities, numerical controls, and declared primary and confirmation
starts. It contains no backend selector.

`RegressionResult` owns:

- resolved problem, Provider, source, and artifact fingerprints;
- ordered fitted parameters with start, final value, unit, movement, bounds,
  scale, and active-bound distance;
- per-row observed/model values, raw/scaled residuals, derivative status, and
  evaluated/skipped/failed accounting;
- Ceres termination and cost diagnostics;
- full and parameter-projected singular values, rank, and conditioning;
- primary/confirmation agreement;
- separate solver, numerical, physical-validity, workflow, scientific, and
  predictive statuses.

Accepted legacy entry points remain stable while their internals migrate to
this engine. Their existing public result types may remain compatibility
projections of the canonical result until an explicitly admitted breaking
change; they must not retain separate solvers or scientific logic. Each
completed migration deletes its superseded paper-specific native cost path.
New paper-named fit functions are forbidden when the paper can be represented
as data plus a `RegressionProblem`.

## Identifiability and numerical preflight

Every observation capability declares its residual and nuisance-variable
dimensions. Before Ceres, Regression verifies that every active parameter has
at least one exact finite derivative path to an observation. A missing or
structurally zero column is a fail-closed input error.

Rank is computed from the scaled Jacobian with

```text
tau = 100 * epsilon_binary64 * max(M, N) * sigma_max
rank = count(sigma_i > tau)
```

For direct observations, the parameter Jacobian is tested directly. For a
lifted formulation with nuisance columns `J_u`, the parameter rank uses
`J_p_projected = (I - Q_u Q_u^T) J_p`, where `Q_u` spans `J_u`. The full lifted
Jacobian rank is also reported. Active-bound coordinates are reported
separately and are not hidden as identified interior parameters.

Primary and confirmation starts receive the same deterministic derivative
preflight. Initial numerical rank deficiency is reported but need not prevent a
nonlinear solve unless all declared starts contain a missing/zero parameter
path. A completed fit cannot be numerically or scientifically valid unless the
final projected parameter rank equals the number of independent active
coordinates. The numerical contract must declare a maximum accepted condition
number and its rationale; the engine supplies no universal scientific cutoff.

## Acceptance and failure rules

A fit is not scientifically successful merely because Ceres terminates.
Results separately report:

- solver convergence;
- numerical convergence and confirmation-start agreement;
- complete exact Jacobian columns;
- local parameter rank and conditioning;
- active bounds and parameter movement;
- observation-domain, phase, topology, and source validity;
- predictive status based only on a predeclared untouched partition and
  approved criteria.

Rank deficiency, unsupported derivatives, Provider fingerprint mismatch,
incomplete row evaluation, invalid physical state, or confirmation failure
returns a diagnostic result or a fail-closed input error according to whether
Ceres began. No runtime automatically persists fitted parameters into a
Provider catalog.

## Performance contract

Broad input support must not move the Ceres evaluation loop into Python.
`RegressionProblem` is validated and serialized once; C++ owns residual
evaluation, Jacobian assembly, Ceres, and diagnostics. Provider calls are
native, model-bound, and batched by compatible topology/observation contract.

The implementation shall:

- tape or prepare a fixed Provider topology once per active problem and reuse
  it when the Provider capability permits;
- reuse row buffers and avoid per-iteration Python objects;
- preserve the exact sparse parameter-to-observation dependency structure;
- choose a deterministic dense or sparse Ceres policy internally from the
  problem shape, with no public backend selector;
- propagate Provider evaluation budgets, deadlines, cancellation, and failed
  row identities;
- retain a benchmark comparing a migrated reference campaign with its legacy
  native path before deleting that path.

Performance evidence cannot relax derivative, rank, or physical-validity
gates.

The installed pair-family campaigns evaluate all 17 audited May 2015
methane/ethane rows as training data. Each `68 x 35` solve converges with full
rank 35 and projected parameter rank 1. The independent fitted values are
non-bound `k_ij = -0.00843032298906253` and
`l_ij = -0.002774426668544412`, each confirmed from two perturbed starts.
This is in-sample reproduction evidence only; the retained pressure-closure
result remains negative under its historical gate and predictive status remains
`NOT_ADJUDICATED_NO_APPROVED_HELD_OUT_CUTOFF`.

The scalar pure implementation replays the four accepted methane training
rows independently for each family. Every `16 x 9` solve has full rank 9 and
projected parameter rank 1, with a non-bound result and confirmation-start
agreement. The local anchors are `m = 1.0001569260577763`,
`sigma = 3.7063548743836034 angstrom`, and
`epsilon/k = 150.00325287725062 K`. They are deterministic in-sample
implementation evidence, not a replacement for the accepted joint fit,
prediction evidence, or Provider-catalog authority.

## Implementation sequence

1. **Complete.** Introduce the typed capability, parameter, observation,
   problem, and result contracts without adding a second native module or
   target.
2. **Complete for fixed-composition VLE.** Generalize neutral,
   non-associating binary `k_ij` from paper-named input ownership to
   caller-supplied component pairs and rows, matching the advertised
   phase-block domain and retaining reference campaigns as evidence.
3. **Complete for fixed-composition VLE.** Add Provider active `l_ij` support
   and admit it through the same phase observation contract and Ceres owner.
4. **Complete for independent scalar pure saturation.** Admit `m`, `sigma`,
   and `epsilon/k` component coordinates through the same contract, result,
   Ceres owner, and native target. The accepted joint methane/ethane workflow
   remains unchanged. Born diameter and solvation factor still require
   migration onto the shared direct-observation path.
5. Add association energy/volume and an explicit `k_hb_ij` combining-rule
   coordinate only after Provider supplies exact topology-bound derivatives
   and source-backed datasets establish rank.
6. Add polar, dielectric, and temperature-dependent families as separately
   evidenced capabilities rather than new engines. Equilibrium-coupled
   observations remain blocked by the dependency doctrine and are not part of
   this sequence.

Each step must be independently fit-ready and reference-validated. The final
goal is broad parameter-family coverage; sequencing controls scientific risk
and provider prerequisites, not the component domain of completed families.

## Explicit exclusions

- arbitrary Python or C++ residual plugins;
- mutable runtime registration;
- copied Provider or Equilibrium equations;
- numerical production derivatives;
- hidden bounds, scales, starts, units, chemistry, or source defaults;
- automatic experimental-uncertainty interpretation;
- parameter-catalog persistence;
- claims of global uniqueness, uncertainty, predictive validity, or
  downstream readiness without their separate evidence.
