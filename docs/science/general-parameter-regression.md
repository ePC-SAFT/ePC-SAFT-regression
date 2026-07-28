# General ePC-SAFT Parameter Regression

Status: user-approved architecture; runtime implementation pending

Date: 2026-07-27

## Objective

Regression shall support source-bound fitting of any admitted continuous
ePC-SAFT parameter for any Provider-supported component set when the caller
supplies an informative dataset and a complete numerical contract.

Published methane, ethane, Figiel, and later literature workflows are
reference campaigns. They prove parameter-family behavior but do not restrict
the runtime to those components, mixtures, or papers.

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
registry. The supported vocabulary is a finite typed schema negotiated with
the installed Provider artifact.

## Capability meanings

The package shall distinguish these states for every parameter family:

1. `REPRESENTED`: Provider can resolve and evaluate a model containing the
   parameter.
2. `DERIVATIVE_READY`: an installed Provider artifact exposes the exact active
   derivative required by at least one observation contract.
3. `FIT_READY`: Regression can construct and solve the corresponding
   source-bound problem for caller-supplied components and rows.
4. `REFERENCE_VALIDATED`: at least one independent source-backed installed-
   artifact campaign has passed its declared gates.
5. `ADMITTED`: an external authority decision, when required, has admitted the
   exact capability. A fitted value is never Provider-catalog authority merely
   because a fit converged.

An example dataset may advance validation evidence without narrowing the
component domain of a `FIT_READY` family.

## Typed parameter identity

An active coordinate is identified by its physical owner, not by a paper or
array position:

```text
ParameterCoordinate
  family
  component_id?                         # single-component family
  component_pair?                      # symmetric pair family
  association_endpoint_pair?           # (component, site) x 2
  correlation_term?                    # named coefficient when applicable
  unit
  start
  lower_bound
  upper_bound
  scale
```

The family and identity fields must match an installed Provider capability.
Units, bounds, scale, and start are mandatory source or workflow inputs; the
runtime supplies no chemistry defaults. Multiple observations may share one
coordinate, and multiple coordinates may be fitted together.

Discrete model choices, charge number, component identity, site topology, and
electrolyte formulation are specifications, not continuous fit coordinates.

## Parameter and data map

The data column names below describe observation families, not guaranteed
identifiability. Every submitted problem must still pass numerical sensitivity
and rank gates.

| Parameter family | Identity | Candidate informative observations | Current state |
|---|---|---|---|
| `segment_count` (`m`) | component | pure saturation pressure/density, PVT density, phase-equilibrium and caloric observations | fit-ready only through the current pure-saturation surface; general request pending |
| `segment_diameter` (`sigma`) | component or declared correlation coefficient | liquid density/PVT, saturation density, phase equilibrium | fit-ready only through the current pure-saturation surface; general request pending |
| `dispersion_energy_over_k` (`epsilon/k`) | component | vapor pressure, phase equilibrium, caloric and PVT observations | fit-ready only through the current pure-saturation surface; general request pending |
| `k_ij` | unordered component pair | VLE/LLE, solubility, activity/fugacity coefficients, MIAC when the electrolyte formulation makes that pair active | bounded Figiel and methane/ethane evidence exists; mixture-independent request pending |
| `l_ij` | unordered component pair | density, excess volume, and phase-equilibrium data sensitive to cross-size mixing | represented by Provider; exact active derivative and Regression surface pending |
| `association_energy_over_k` | association endpoint pair | hydrogen-bond-sensitive VLE/LLE, solvation, enthalpy, and related phase-property observations | represented by Provider; exact active derivative and Regression surface pending |
| `association_volume` | association endpoint pair | hydrogen-bond-sensitive VLE/LLE, density, solvation, and caloric observations | represented by Provider; exact active derivative and Regression surface pending |
| `k_hb_ij` | source-defined cross-association combining-rule coordinate | the same cross-association-sensitive observations, with enough composition/temperature variation to separate it from pure association parameters | not a current Provider record family; requires an explicit active combining-rule coordinate or a declared transform to resolved cross-association energy |
| `born_diameter` | ion component | source-defined single-ion solvation Gibbs energy or another explicitly admitted ionic observable | five-ion fit demonstrated; general request pending |
| `solvation_factor` | applicable component | MIAC, osmotic/activity, or solvation observations under a declared reference sequence | one NaBr water-factor fit demonstrated; general request pending |
| dielectric and ion-suppression coefficients | component/model/correlation term | dielectric, MIAC, osmotic/activity, solvation, or phase observations over a rank-sufficient state range | represented families vary; exact active derivatives and Regression surface pending |
| polar coefficients | component/correlation term | polar-mixture VLE, PVT, caloric, dielectric, or other source-defined polar observables | represented families vary; exact active derivatives and Regression surface pending |
| temperature-dependent coefficients | named correlation term | observations spanning enough temperatures to identify every active coefficient | represented correlations vary; exact active derivatives and Regression surface pending |

The table describes physically relevant evidence, not an automatic weighting
policy. Experimental uncertainty may define a weight only when the source
actually reports it and the request declares that interpretation.

### `k_hb_ij` rule

`k_hb_ij` must not be inferred from a resolved cross-association energy and
silently fitted under another name. A fit may use it only when Provider exposes
the source-defined combining rule, its component/site identity, and its exact
total derivative. Otherwise the caller may fit the resolved
`association_energy_over_k` coordinate, and the result must retain that
different meaning.

## Source-bound observations

Every observation row contains:

```text
row_id
observation_type
component_ids and exact order
T and P specifications
composition and basis
phase or reference-state convention
observed value and unit
residual scale and its rationale
partition
source_id, locator, transformation, and artifact hash
```

Initial typed observation contracts are:

- pure saturation coexistence;
- single-phase scalar property or density;
- fixed-composition phase equilibrium with observed phase compositions and
  lifted phase state variables;
- activity, fugacity, osmotic, or mean-ionic-activity coefficient;
- solvation Gibbs energy.

Positive equality, linear aggregate, and one-sided censored observations may
be added as typed contracts when a real source requires them. Equilibrium-
eliminated observations require an installed Equilibrium value and exact
implicit sensitivity contract; Regression must not run a second equilibrium
implementation or finite-difference a black-box solve.

The loader rejects missing units, ambiguous composition bases, duplicate row
identities, unsupported phases, nonfinite values, missing source identities,
and parameter/observation combinations that the installed artifact does not
advertise.

## Exact derivative contract

The installed Provider artifact shall expose a finite capability description
and model-bound evaluators. Each capability identifies:

- parameter family and identity shape;
- accepted observation or phase-block contract;
- active coordinate order and units;
- returned value units;
- exact derivative order;
- fixed topology and domain limits;
- Provider parameter and artifact fingerprints.

Direct observable residuals consume Provider values and first total
derivatives:

```text
r_q = (y_q(p) - y_q,observed) / s_q
J_qj = (dy_q/dp_j) * parameter_scale_j / s_q
```

Lifted phase-equilibrium residuals consume Provider phase-potential values,
gradients, and Hessians. The Hessian supplies exact derivatives of pressure
and chemical-potential residuals with respect to lifted state variables and
active parameters; no third derivative is required for the Ceres Jacobian.

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
- separate solver, numerical, workflow, scientific, and predictive statuses.

Accepted legacy entry points remain stable while their internals migrate to
this engine. New paper-named fit functions are forbidden when the paper can be
represented as data plus a `RegressionProblem`.

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

## Implementation sequence

1. Introduce the typed capability, parameter, observation, problem, and result
   contracts without adding a second native module or target.
2. Generalize `k_ij` from paper-named input ownership to caller-supplied
   component pairs and rows, retaining the existing reference campaigns as
   installed-artifact evidence.
3. Add Provider active `l_ij` support and admit it through the same phase and
   property observation contracts.
4. Migrate pure `m`, `sigma`, and `epsilon/k`, Born diameter, and solvation
   factor workflows onto the shared contracts while preserving accepted
   methane and ethane numerical behavior.
5. Add association energy/volume and an explicit `k_hb_ij` combining-rule
   coordinate only after Provider supplies exact topology-bound derivatives
   and source-backed datasets establish rank.
6. Add polar, dielectric, temperature-dependent, and equilibrium-coupled
   families as separately evidenced capabilities rather than new engines.

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
