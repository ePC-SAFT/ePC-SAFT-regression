# General ePC-SAFT Parameter Regression

Status: ten nonpolar/electrolyte parameter families fit-ready for their typed
observation contracts; pure 2B association derivative/surface-ready but
source-series blocked; polar regression excluded

Date: 2026-07-28

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
| `rueben_dipole_scaling`, `rueben_polarizability_scaling`, `rueben_correlation_integral_parameter` | Provider-represented polar-model inputs explicitly excluded from the Regression roadmap by user decision; no Regression capability is planned |
| `k_ij`, `l_ij` | continuous unordered component-pair coordinates |
| `association_energy_over_k`, `association_volume` | continuous model-bound coordinates only for an advertised pure 2B model with one symmetric association pair; general association-endpoint identity remains pending |
| `dielectric_ion_suppression_coefficient`, `ionic_region_relative_permittivity` | continuous model coordinates when exact derivatives are advertised |
| correlation terms for `segment_diameter`, `relative_permittivity`, and `solvation_factor` | individually named coefficient coordinates; no whole-correlation opaque fit |

New Provider record families require an explicit schema revision in both
packages. Regression must never infer fit support merely because a record can
be resolved.

## Definitive parameter-family and data-requirement matrix

`FIT_READY` below means that caller-supplied rows can execute today when the
installed Provider advertises the named exact capability. It does not mean
that every dataset is informative, that a reference value will be reproduced,
or that the result is predictive. Every problem still needs declared units,
source identity, partitions, residual scales, parameter bounds/scales/starts,
complete row accounting, full projected parameter rank, acceptable
conditioning, non-bound diagnostics, and confirmation-start agreement.

| Parameter family | Identity | Minimum source-bound data contract | Exact executable contract | Current state and MEA role |
|---|---|---|---|---|
| `segment_count` (`m`) | component | Multiple pure rows with `T`, observed saturation pressure and liquid density, molar mass, phase-volume bounds/starts, partitions, and scales | Provider `(n,V,m)` value/gradient/Hessian; Regression lifted liquid/vapor volumes with pressure, chemical-potential, and density residuals | `FIT_READY`, one family at a time. Standalone recovery is optional for MEA unless `m` is selected in MEA's application-owned parameter block. |
| `segment_diameter` (`sigma`) | component | Same pure saturation pressure/liquid-density contract, over a range that gives nonzero independent sensitivity | Provider `(n,V,sigma)` value/gradient/Hessian; same lifted Ceres owner | `FIT_READY` for a constant coordinate. Named temperature-correlation coefficients are separately `NOT_READY`. Optional for MEA unless selected. |
| `dispersion_energy_over_k` (`epsilon/k`) | component | Same pure saturation pressure/liquid-density contract with enough vapor-pressure sensitivity | Provider `(n,V,epsilon/k)` value/gradient/Hessian; same lifted Ceres owner | `FIT_READY`, one family at a time. Optional for MEA unless selected. |
| `k_ij` | unordered component pair | One supported domain: fixed measured `T,P,x,y` VLE rows; source-bound MIAC rows with formula molality; or single-ion solvation-Gibbs targets. Rows must vary enough to identify the selected pair. | Exact Provider Hessian for lifted neutral VLE or exact first derivative for the admitted aqueous/solvation direct observable; one typed pair coordinate | `FIT_READY` in the advertised neutral-VLE, aqueous-MIAC, and organic-ion-solvation domains. MEA needs only the exact pair derivatives for pairs actually selected in its 12-parameter block. |
| `l_ij` | unordered component pair | Fixed measured `T,P,x,y` VLE rows sensitive to cross diameter, with explicit source sign convention | Exact Provider `(n1,n2,V,l_ij)` Hessian for the currently admitted neutral nonassociating binary domain | `FIT_READY` only for fixed-composition neutral VLE. Density, excess-volume, associating, and electrolyte observation domains are `NOT_READY`; standalone expansion is not an MEA prerequisite unless selected. |
| `born_diameter` | ion component | One or more source-defined single-ion solvation-Gibbs targets with exact x-process convention, state, component order, and numerical scale | Exact Provider solvation-Gibbs value/first derivative for the active ion; direct-observable Ceres row | `FIT_READY`; five Figiel ions are reference evidence. It is parallel parameter groundwork, not coupled-MEA readiness. |
| `solvation_factor` | component | Source-bound MIAC rows at declared `T,P` and formula-unit molality, with solvent/ion identities and scale | Exact Provider MIAC value/first derivative for one active factor; direct-observable Ceres rows | `FIT_READY` for the advertised constant factor. A temperature correlation is separately `NOT_READY`; not on the MEA critical path unless selected. |
| `relative_permittivity` | solvent component | Source-bound single-ion solvation-Gibbs targets with fixed other model inputs, exact solvent identity, x-process convention, state, and scale | Exact Provider solvation-Gibbs value/first derivative for the active solvent permittivity; direct-observable Ceres rows | `FIT_READY`. Corrected Validation subject `e4cb7af` gives five independent rank-1 water fits returning `78.0899937514462` through `78.08999375166104` versus fixed `78.09`. This is implementation evidence, not a paper-fitted target or MEA prerequisite. |
| `dielectric_ion_suppression_coefficient` | model | Salt-free-normalized relative-permittivity observations spanning enough total-ion mole fraction to identify one coefficient | Exact Provider relative-permittivity ratio and first derivative; direct-observable Ceres rows | `FIT_READY`; 36 digitized Figiel water/methanol rows are reference evidence. Optional standalone recovery. |
| `ionic_region_relative_permittivity` | model | Source-defined single-ion solvation-Gibbs targets with fixed Born diameters and other inputs | Exact Provider SSM+DS solvation-Gibbs value/first derivative; direct-observable Ceres rows | `FIT_READY`; five independent Figiel fits recover fixed `8` within `2.2e-9`. Parallel groundwork, not coupled-MEA readiness. |
| `association_energy_over_k` | association endpoints; currently one model-bound symmetric pure 2B pair | Simultaneous multi-temperature vapor-pressure and liquid-density observations; one density point is insufficient. The first campaign uses the frozen Baygi--Pahlavanzadeh MEA target correlations. | Current scalar `(n,V,p)` Hessian is mechanics evidence only. The selected route requires one joint Provider Hessian and the existing lifted pure-saturation owner over `(n,V,m,sigma,epsilon/k,epsilon_assoc/k,kappa_assoc)`. | `BAYGI_MEA_2B_DESIGN_READY_PROVIDER_JOINT_CALLBACK_REQUIRED`. Gross--Sadowski and Baygi--Pahlavanzadeh did not recover this coordinate independently. |
| `association_volume` | association endpoints; currently one model-bound symmetric pure 2B pair | The same simultaneous vapor-pressure/liquid-density series is required to distinguish this coordinate from association energy and the three ordinary pure parameters | Same selected joint five-parameter Provider and Regression contract; no second association-only engine | `BAYGI_MEA_2B_DESIGN_READY_PROVIDER_JOINT_CALLBACK_REQUIRED`. |
| `k_hb_ij` | source-defined cross-association combining-rule coordinate | A source-defined combining rule and sign convention plus composition/temperature-varying association-sensitive observations that separate cross from pure association | New Provider record/transform identity, versioned fingerprint, exact chain-rule derivative, and an association-endpoint Regression identity are required | `NOT_READY`. Ascani's fixed `0.026` is provenance, not a recovery dataset. Do not alias it to resolved association energy or invent zero defaults. Not the next MEA investment. |
| `schreckenberg_dielectric_volume`, `schreckenberg_dielectric_temperature` | component/correlation coordinates | Multi-temperature electrolyte relative-permittivity or other source observables that independently identify the selected coefficient | New exact Provider active-parameter callback and Regression correlation identity/surface | `REPRESENTED_NOT_DERIVATIVE_READY`; no retained rank-sufficient recovery series. Optional, not on the MEA critical path. |
| `zuber_ion_suppression_coefficient` | ion component | Relative-permittivity, MIAC, osmotic/activity, or solvation observations spanning ion fraction/molality and preferably temperature | New exact Provider active-parameter callback and matching typed Regression observation contract | `REPRESENTED_NOT_DERIVATIVE_READY`; the retained one-row osmotic oracle is insufficient. Optional. |
| Named correlation coefficients for `segment_diameter`, `relative_permittivity`, or `solvation_factor` | component, correlation family, coefficient kind, and term index | Multi-temperature rows covering a range sufficient to distinguish constant, amplitude, and exponent coefficients | New closed correlation-coordinate descriptor, exact Provider derivative including the correlation chain rule, and Regression identity | `REPRESENTED_NOT_DERIVATIVE_READY`; current constant-parameter evidence does not validate coefficient recovery. |
| Temperature-dependent `k_ij(T)` coefficients | pair plus named correlation coefficient | Multi-temperature mixture observations and the exact source correlation form, reference temperature, units, and sign convention | New Provider pair-correlation record/capability and exact derivative; Regression correlation identity | `NOT_REPRESENTED_NOT_READY`. Ascani provides a source law but the retained bounded case is only 298.15 K. |
| Polar parameter families | source-specific | none in this roadmap | none in Regression | `EXCLUDED_FROM_REGRESSION_ROADMAP_BY_USER_DECISION`. Provider may still represent and evaluate polar physics for other owners. |

`molar_mass`, `charge_number`, component identity, association-site topology,
association scheme, and electrolyte formulation remain fixed measured or
discrete model inputs. They are not continuous Regression coordinates.
Experimental uncertainty may define a residual weight only when the source
reports it and the submitted contract declares that interpretation.

## Downstream decision: MEA before optional family recovery

Standalone parameter recovery and coupled reactive regression are different
deliverables. The Figiel/Born, pure-association, dielectric, and correlation
campaigns can demonstrate individual Provider/Regression contracts, but they
do not establish the coupled state sensitivities needed by
MEA-Thermodynamics.

MEA's critical Regression path is:

1. Migration's D-026 installed two-liquid Stage-II/III gate passes. This
   remains the required public prerequisite; the MEA sequence must not infer
   admission from private reacting-phase foundations or skip the installed
   two-liquid evidence.
2. Provider exposes exact parameter partials for only the parameter identities
   selected by the MEA application. Regression must not infer a broader
   chemistry block or persist fitted values to the Provider catalog.
3. Equilibrium returns converged reactive-liquid and reactive-bubble values
   together with exact implicit sensitivities to those parameters. The value
   solve remains outside the tape and the sensitivity contract follows
   `u_z = -H_u^{-1} H_z`, with conditioning and failure diagnostics. Regression
   must not copy reaction/EOS equations, run a second equilibrium formulation,
   or finite-difference a black-box solve.
4. Regression reuses its one native Ceres engine/result/target and adds only
   the schema-driven parameter sharing plus three observation semantics that
   real MEA evidence requires: positive equality, linear aggregate, and
   one-sided censored observations. Every row retains its ID, provenance,
   partition, scale/weight or censor policy, and evaluated/skipped/failed
   accounting with complete Jacobian columns.
5. The first installed falsification slice is the reduced two-row MEA fit: one
   accepted reactive-liquid CO2-loading equality and one accepted
   reactive-bubble heat-of-absorption/enthalpy equality. It must pass exact
   Jacobian, rank/conditioning, bound/KKT, solver, numerical, equilibrium
   physical-validity, and row-accounting gates before broadening.
6. Only then may the application-owned frozen campaign use MEA's 9 species,
   5 reactions, 12-parameter block, starts/bounds/regularization, 147-state
   training partition with 297 observations, and untouched 220-state reserved
   partition with 435 observations. Reserved states are evaluated without
   refitting; application-owned promotion cutoffs remain separate.

Therefore:

`NEXT_ENGINEERING_INVESTMENT_AFTER_D026: EQUILIBRIUM_EXACT_IMPLICIT_PARAMETER_SENSITIVITY_FOR_MEA`

This is chosen over acquiring the missing pure/cross-association recovery
series. Those sources remain useful for optional standalone validation, but
they do not unblock the first coupled MEA tracer. The immediate implementation
must still wait for D-026, an exact installed Equilibrium value/sensitivity
artifact, and the matching Provider parameter-partial contract; this document
does not authorize speculative Regression runtime.

### `k_hb_ij` rule

`k_hb_ij` must not be inferred from a resolved cross-association energy and
silently fitted under another name. A fit may use it only when Provider exposes
the source-defined combining rule, its component/site identity, a versioned
transform identifier and fingerprint, and its exact total derivative including
the transform chain rule. Otherwise the caller may fit the resolved
`association_energy_over_k` coordinate, and the result must retain that
different meaning.

### Pure-association evidence and selected joint route

Provider commits `c9ada20` and `a4d8a0e` add no new callback or solver. They
extend the existing scalar pure phase callback only when the installed model
has one neutral component, exactly two association sites, and one symmetric
active pair. The descriptor identity is model-bound and its topology
fingerprint is mandatory. It therefore does not claim arbitrary association
endpoint, mixed association, or `k_hb_ij` support.

For one pure-density observation, Regression owns the lifted coordinate

```text
p = p_origin + p_scale z
V = V_origin exp(u)
```

and the two scaled residuals

```text
r_P   = (P_model(T, n=1, V, p) - P_observed) / s_P
r_rho = (M / V - rho_observed) / s_rho
```

The exact Jacobian consumes the Provider Hessian in `(n,V,p)`:

```text
dr_P/dz   = -R T Phi_Vp p_scale / s_P
dr_P/du   = -R T Phi_VV V / s_P
dr_rho/dz = 0
dr_rho/du = -(M/V) / s_rho
```

Thus one row is a `2 x 2` local problem and must pass full rank 2 plus
projected parameter rank 1. This surface is executable with the direct
Held-2012 ethanol density anchor retained by Validation (CSV SHA-256
`25e3be94ee3cfb5eb13df89827f4368673f87f96d6b0225c12e35f9396b8779c`),
using 1 bar as the declared model approximation to the reported ambient
pressure, but that single point is not the paper's parameter-recovery target.
For both independent parameter families, the two declared one-row starts
return finite full-rank diagnostics but do not reach solver or numerical
convergence. That is retained surface/falsification evidence, not a fitted
parameter result.

Gross and Sadowski (2002, DOI `10.1021/ie010954d`) state that all five pure
parameters were adjusted simultaneously against vapor-pressure and
liquid-density data over the reported temperature range; for ethanol the
range is 230--516 K. The read-only Markdown artifact used to verify that
method has SHA-256
`dc4695f03a2511f0ac416bfb54923ed2b7b7a9ced8240d10b112b42ad977d732`.
Table 1 supplies the fitted coordinates, aggregate AAD values, temperature
ranges, and source-compilation references, but no raw rows, row counts,
objective equation, weights, starts, bounds, optimizer, or stopping criteria.
The detailed evidence audit is
`docs/research/gross-sadowski-2002-association-fitting.md`. Consequently the
two association families are not independently `FIT_READY`.

The selected implementation route is one extension of the existing
pure-saturation Ceres owner, not a new association fitter. If the acquired
source packet contains `N` vapor-pressure/liquid-density observations at
common temperatures, it has five global transformed parameters

```text
(m, sigma, epsilon/k, epsilon_assoc/k, kappa_assoc)
```

and two lifted log volumes per row. The problem therefore has `5 + 2N`
variables and `4N` residuals: liquid pressure closure, vapor pressure closure,
liquid--vapor chemical-potential equality, and observed liquid density. The
exact Jacobian requires one Provider value/gradient/Hessian block in

```text
(n, V, m, sigma, epsilon/k, epsilon_assoc/k, kappa_assoc)
```

for each phase. No third derivatives, density-root callback, copied EOS,
numerical derivative backend, second Ceres engine, or association-only result
family is required. Admission requires full rank `5 + 2N`, projected
parameter rank 5, acceptable conditioning, non-bound diagnostics, and
independent-start confirmation.

Experimental-row campaigns must preserve their source temperatures rather
than interpolate them. A correlation-defined campaign may predeclare a
deterministic evaluation grid, but it must identify those rows as calculated
correlation targets rather than experimental measurements. If a later source
packet contains independent pressure and density grids, the same owner must
retain them as separate typed observations and derive its dimensions from
those rows; it must not manufacture paired points merely to reuse the current
`4N` shape.

This formulation reproduces the papers' verified simultaneous parameter
coupling. Its residual scaling, bounds, starts, and row selection must be
predeclared as reconstruction choices unless an authoritative source for an
author's exact numerical method is acquired. The single Held-2012 density
anchor remains derivative/mechanics evidence and must not be used to tune the
joint fit.

### First campaign: Baygi--Pahlavanzadeh MEA 2B reconstruction

The first joint campaign is
`baygi-pahlavanzadeh-2015-mea-2b-correlation-reconstruction-v1`.
The primary source is Baygi and Pahlavanzadeh, *Chemical Engineering
Research and Design* 93 (2015) 789--799,
DOI `10.1016/j.cherd.2014.07.017`. The inspected PDF has SHA-256
`7e8e77577a34bd9867489faee992dd192e8cbbc728c50a26e8264b0e09192365`.
The component identity is monoethanolamine, component id
`monoethanolamine`, CAS `141-43-5`, formula `C2H7NO`, and molar mass
`0.0610831 kg/mol` from NIST Chemistry WebBook SRD 69. The fixed association
topology is 2B: one donor site and one acceptor site with one symmetric pure
association pair.

Baygi and Pahlavanzadeh Table 1 and Eqs. 9--10 define the calculated targets:

```text
P_sat(T) [Pa]
  = exp(92.624 - 10367/T - 9.4699 ln(T) + 1.9e-18 T^6)

rho_L,sat(T) [mol/L]
  = 1.0011 / 0.22523^[1 + (1 - T/678.2)^0.21515]
```

Here `T` is in kelvin. The density form is the standard DIPPR-105 grouping.
It is also the only grouping consistent with the paper's Figure 3: it gives
approximately `1008.5 kg/m3` at `303.15 K` and `889.6 kg/m3` at `443.15 K`;
the literal alternative placement of the printed superscript gives about
`5 mol/L` and contradicts that figure. The implementation source packet must
record this disambiguation and convert mol/L to kg/m3 using the stated molar
mass.

The paper states a fit range of `303.15--443.15 K` but does not publish
`np` or the evaluation grid. This reconstruction therefore freezes 15
calculated training rows:

```text
T_j = 303.15 K + 10 j K,  j = 0,...,14.
```

All 15 rows are training targets. There is no held-out, stress, predictive, or
experimental-uncertainty claim. A later five-kelvin grid replay may measure
grid sensitivity, but it is validation evidence and cannot silently replace
the canonical grid.

The five global parameter coordinates, in fixed order, are:

```text
(m, sigma [angstrom], epsilon/k [K],
 epsilon_assoc/k [K], kappa_assoc [1]).
```

The campaign bounds, affine scales, primary start, and confirmation start are
reconstruction choices selected before execution:

| Coordinate | Bounds | Scale | Primary start | Confirmation start |
|---|---:|---:|---:|---:|
| `m` | `[0.5, 5.0]` | `0.5` | `2.5` | `3.25` |
| `sigma` | `[2.0, 5.0] angstrom` | `0.5 angstrom` | `3.5` | `3.0` |
| `epsilon/k` | `[50, 400] K` | `50 K` | `225` | `300` |
| `epsilon_assoc/k` | `[250, 5000] K` | `500 K` | `2000` | `3000` |
| `kappa_assoc` | `[0.001, 0.25]` | `0.05` | `0.05` | `0.10` |

These bounds enclose the associating parameter sets tabulated by the inspected
Gross--Sadowski, Baygi--Pahlavanzadeh, Diamantonis--Economou, and Fuenzalida
sources. They are campaign bounds, not universal PC-SAFT limits. Neither start
is the published MEA 2B tuple or a rounded copy of it.

For `N=15`, the existing lifted-volume problem has 35 variables and 60
residuals:

```text
variables = five transformed global parameters
          + (log V_L,j, log V_V,j) for j=1,...,15

r_j = [
  (P_L - P_obs) / P_obs,
  (P_V - P_obs) / P_obs,
  mu_L/RT - mu_V/RT,
  (rho_L - rho_obs) / rho_obs
]
```

The four row residuals retain the existing equal `0.25` weights. This is the
smallest extension of the accepted Ceres workflow, but it is not Baygi and
Pahlavanzadeh's exact Eq. 8 objective: their printed objective is an
unweighted sum of absolute relative errors in the final saturation pressure
and liquid density, and their optimizer is not reported. The result must
therefore calculate and report the two paper-style AAD values after solving
the equilibrium reporting problem, while labeling the fitted result a
correlation reconstruction rather than an exact author-run replay.

The Provider prerequisite is one model-bound callback over the exact
coordinate order

```text
(n, V, m, sigma, epsilon/k, epsilon_assoc/k, kappa_assoc)
```

returning `Phi`, its gradient of length 7, and its symmetric Hessian of shape
`7 x 7`, plus pressure, `mu/RT`, stability diagnostics, parameter fingerprint,
and topology fingerprint. Regression consumes the Hessian entries
`Phi_VV`, `Phi_Vp`, `Phi_nV`, and `Phi_np` for every active parameter `p`.
The exact residual Jacobian follows by the existing chain rule:

```text
dP/dp       = -R T Phi_Vp
dP/d(log V) = -R T Phi_VV V
d(mu/RT)/dp = Phi_np
d(mu/RT)/d(log V) = Phi_nV V
d(rho)/d(log V) = -rho
```

No third derivatives are required. The callback must work with a
`ParameterBundle.from_records` MEA 2B input; the campaign must not persist the
start or fitted values into the Provider catalog.

The same native target, Ceres engine, Python workflow, and
`PureSaturationFitResult` remain the sole owners. Their parameter tuples and
row counts become data-sized rather than adding a second association fitter.
Methane, ethane, and propane records and numerical behavior remain unchanged.

Acceptance of the executable reconstruction requires:

- Provider value/gradient/Hessian agreement with an independent directional
  finite-difference oracle at representative liquid and vapor states;
- finite exact `60 x 35` residual Jacobian with directional agreement;
- full Jacobian rank 35;
- projected parameter rank 5 after eliminating the lifted-volume column
  space, using `(I - J_V J_V^+) J_p`, with retained singular values and
  condition number;
- finite interior fitted parameters and mechanically stable, correctly
  ordered liquid and vapor states;
- converged primary and confirmation starts with retained scaled parameter and
  cost deltas;
- separate solver, numerical, physical/workflow, scientific-comparison, and
  predictive statuses; and
- reported pressure and density AADs plus absolute and relative differences
  from Baygi and Pahlavanzadeh's descriptive 2B tuple
  `(3.0353, 3.0435 angstrom, 277.174 K, 2586.3 K, 0.037470)` and AADs
  `(0.62%, 0.12%)`.

No parameter-distance or AAD acceptance cutoff is invented. A converged,
full-rank reconstruction can establish that the package supports a joint
five-parameter pure 2B fit and can describe its agreement with the published
case. It cannot establish exact reproduction of the authors' undisclosed
optimization, predictive validity, global uniqueness, uncertainty, or
Provider-catalog authority.

Rejected first-slice alternatives are:

1. An exact Eq. 8 author-run replay. It is not currently definable because the
   paper omits its grid, optimizer, starts, bounds, and tolerances; implementing
   a new nonsmooth inner-equilibrium objective would broaden the solver surface
   without recovering those facts.
2. The Albers--Sadowski five-point PCP-SAFT minimum-data design. It is useful
   identifiability evidence, but it is a different model and target selection,
   not the Baygi MEA campaign.
3. An association-only two-parameter fit. The inspected sources fit all five
   coordinates jointly, and the current one-row scalar surface cannot identify
   the source problem.

Electrolyte MEA, aqueous MEA, reaction equilibria, binary `k_ij`,
cross-association, and `k_hb_ij` remain later, separately source-bound
families. This campaign supplies a fitted pure MEA 2B parameter artifact as
optional groundwork; it does not admit or design those later residuals.

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

For aqueous mean-ionic-activity fits, one scalar `k_ij` is active per problem.
Each row carries the ordered Provider model identity and the explicit fixed
`(k_water,cation, k_water,anion, k_cation,anion)` context. The active pair
replaces exactly one entry before the installed Provider batch callback is
evaluated. The residual and exact Jacobian are

```text
r_i = (1 - gamma_model,i / gamma_observed,i) / s_i
dr_i/dz = -(gamma_model,i / gamma_observed,i)
           * d ln(gamma_model,i)/d k_ij * (d k_ij/dz) / s_i.
```

The fixed context is a required workflow input, not a catalog default inferred
by Regression. Water--cation, water--anion, and cation--anion are separate
closed Provider capabilities; fitting one does not silently activate the other
two.

For an advertised organic ion-solvation capability, every row carries the
ordered `(solvent, cation, anion)` model, active ion, active unordered pair,
and explicit fixed
`(k_solvent,cation, k_solvent,anion, k_cation,anion)` context. The active pair
replaces exactly one entry. With solver coordinate `z`, physical pair
parameter `k = k_origin + k_scale z`, and solvation-energy scale `s_G`,

```text
r_G = (G_solv(k) - G_solv,observed) / s_G
dr_G/dz = (d G_solv/d k_ij) * k_scale / s_G.
```

Provider owns the pure-solvent infinite-dilution reference sequence and only
returns success after its value and derivative limits converge. Regression
does not copy the reference sequence, Born term, or EOS. Solvent--cation,
solvent--anion, and cation--anion coordinates are separate closed capability
IDs. The active component must be the cation or anion, never the solvent.

For the model-level Figiel dielectric ion-suppression fit, every row supplies
the total ion mole fraction and the observed dimensionless ratio to the
salt-free solvent permittivity. The residual and exact Jacobian are

```text
r_i = (epsilon_model,i / epsilon_saltfree,i - y_i) / s_i
dr_i/dz = (d epsilon_model,i / d alpha)
           * (d alpha/dz) / (epsilon_saltfree,i * s_i).
```

Provider supplies both permittivities and the exact first derivative;
Regression does not copy the dielectric equation. Solvent identity remains
source provenance even where normalization makes the scalar model response
independent of the salt-free solvent value.

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
- salt-free-normalized relative permittivity.

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
`(n,V,active_parameter)` value/gradient/Hessian callback. Provider
`2d1816cf376294156684fee85611a93fc41d0970` adds complete typed descriptors
and topology fingerprints for the existing active-Born and aqueous-
solvation-factor callbacks. Regression selects the requested known descriptor
from that finite Provider capability set, reports unknown
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

The direct-observable descriptors have derivative order one because Provider
returns the observable and its exact total derivative with respect to the sole
active physical parameter. Regression's typed direct-observation contracts
make those exact domains `FIT_READY`; this does not admit fitted values to the
Provider catalog or generalize their fixed reference paths.

Direct observable residuals consume Provider values and first total
derivatives. For a solvation-Gibbs target, solver coordinate `z`, physical
parameter `p = p_origin + p_scale z`, and residual scale `s_G`:

```text
r_G = (G_solv(p) - G_solv,observed) / s_G
J_G = (dG_solv/dp) * p_scale / s_G
```

The mean-ionic-activity contract preserves the source-frozen relative
residual rather than substituting a log residual:

```text
gamma_model = exp(log_gamma_model)
r_gamma = (1 - gamma_model/gamma_observed) / s_relative
J_gamma = -(gamma_model/gamma_observed)
          * d(log_gamma_model)/dp * p_scale / s_relative
```

The implemented direct path accepts exactly one active parameter per problem.
Parameter sharing is therefore explicit across rows, while the five Figiel
Born targets are five independent one-parameter problems. Organic ion-solvent
pair fits are likewise independent; other pair values remain explicit row
inputs.

The reference organic-ion campaign uses four constructed nearest-pure
endpoints retained by Validation: K+/methanol, Br-/methanol, Na+/ethanol, and
Cl-/ethanol. Three digitized compositions are only near unity. Therefore this
campaign proves exact-derivative fit execution and provides a descriptive
comparison with rounded Table 5 parameters; it is not exact reconstruction of
the paper's pure-organic fitting data, and the remaining mixed-solvent rows
are not duplicated as training residuals.

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
   remains unchanged.
5. **Complete for the advertised direct-observable domains.** Admit one Born
   diameter from one or more solvation-Gibbs targets, or one solvation factor
   from one or more mean-ionic-activity targets. Each problem remains one
   parameter at a time, consumes the Provider value/first derivative, and
   uses no lifted phase variables.
6. **Complete for the advertised model-level dielectric domain.** Admit one
   ion-suppression coefficient from normalized relative-permittivity rows.
   The 36-row Figiel water/methanol reference fit is a `36 x 1`, rank-1,
   non-bound training solve and recovers `7.067350349980952` versus the
   paper's descriptive `7.01`; the digitized rows define no scientific or
   predictive cutoff.
7. **Complete for model-level ionic-region relative permittivity.** Admit one
   positive SSM+DS coordinate from source-bound solvation-Gibbs targets on an
   installed model-bound callback. Validation subject `81e9a3f` fits all five
   admitted Table S5 ion models independently; every `1 x 1` problem is rank
   1, non-bound, and confirmed from starts 4 and 12, recovering the fixed
   input `8` within `2.2e-9`. This is in-sample mechanics evidence, not
   recovery of a parameter that Figiel fitted or a shared multi-model value.
8. **Complete for component solvent relative permittivity.** Admit one
   positive solvent coordinate from source-bound single-ion solvation-Gibbs
   targets. Provider subject `7cadaad`, Regression subject `177890a`, and
   corrected Validation subject `e4cb7af` bind the exact installed-artifact
   path. Five
   independent `1 x 1` fits are rank 1, non-bound, and confirmed from starts
   50 and 110. Recovery of fixed `78.09` is mechanics evidence, not a
   paper-fitted or predictive result.
9. **Partial for pure 2B association.** Exact Provider Hessians and the
   Regression pure-density surface exist for the sole symmetric association
   pair. Recovery remains source-blocked until the original simultaneous
   vapor-pressure/liquid-density rows and objective are retained. General
   association endpoints and an explicit `k_hb_ij` combining-rule coordinate
   remain separate future capabilities.
10. **Deferred optional standalone recovery.** Remaining nonpolar dielectric,
    ion-suppression, association, and temperature-correlation families require
    their own source series plus exact Provider derivative seams. They do not
    justify speculative Regression implementation. Polar families are
    excluded from this roadmap.
11. **Selected next investment: coupled MEA prerequisite.** After the required
    D-026 installed two-liquid Stage-II/III gate, obtain an exact installed
    Equilibrium reactive-state value/implicit-parameter-sensitivity contract,
    paired with exact Provider parameter partials, before adding the reduced
    two-row mixed-observable Ceres tracer.

Every standalone family admitted in steps 1--10 must be independently
fit-ready and reference-validated. Step 11 is a separate coupled capability:
it does not promote optional standalone families and cannot use their evidence
as a substitute for exact Equilibrium sensitivities.

## Explicit exclusions

- arbitrary Python or C++ residual plugins;
- mutable runtime registration;
- copied Provider or Equilibrium equations;
- numerical production derivatives;
- hidden bounds, scales, starts, units, chemistry, or source defaults;
- automatic experimental-uncertainty interpretation;
- parameter-catalog persistence;
- polar parameter regression;
- claims of global uniqueness, uncertainty, predictive validity, or
  downstream readiness without their separate evidence.
