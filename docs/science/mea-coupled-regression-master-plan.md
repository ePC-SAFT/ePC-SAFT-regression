# MEA Coupled Regression Master Plan

Status: source-bound coordination contract; runtime not ready.

This document is the single Regression-owned plan for fitting ePC-SAFT
parameters to aqueous monoethanolamine observations. It coordinates the
application contract owned by MEA-Thermodynamics, thermodynamic state and
sensitivity results owned by Equilibrium, EOS values and exact parameter
partials owned by Provider, and the Ceres workflow owned by Regression.

It does not admit a parameter set, authorize a fit, persist values to the
Provider catalog, or make a predictive claim.

## 1. Governing decision

The first coupled tracer shall use:

1. one admitted CO2 partial-pressure observation; and
2. one eligible liquid-speciation equality or linear aggregate.

CO2 loading is a fixed state input, not the measured response. Heat of
absorption is excluded from the first tracer because no source-complete heat
row or observation definition is admitted.

Under the current MEA application contract, the pressure prediction is

```text
p_CO2 = y_CO2 * P_bubble
```

from a certified reactive-bubble calculation with Provider fugacity equality.
The bubble result must additionally certify a zero minimized Provider
phase-potential/TPD condition, vapor incidence and normalization, transfer
equilibrium, positivity, packing/domain validity, and an independently
recomputed final state. Multiple roots, finite-search completion, and
`globality = not_guaranteed` remain separate diagnostics. No ideal-vapor
shortcut is frozen. Replacing this contract with
`p_CO2 = f_CO2_liquid` would require a new MEA-owned source and domain
justification; Regression shall not infer that approximation.

The old proposed tracer—CO2-loading equality plus heat equality—is rejected as
unexecutable from the admitted source packet. The old historical 12-parameter
block and 435 reserved-observation count are also stale.

## 2. Ownership

| Owner | Owns | Must not own |
|---|---|---|
| MEA-Thermodynamics | species and reactions; equilibrium-constant sources; standard-state convention; feed/state construction; source rows; target eligibility; measurement mapping; parameter selection; bounds, scales, starts, sharing, regularization, partitions, and scientific gates | Provider equations, Equilibrium algorithms, Ceres, or fitted-value persistence |
| Provider | typed parameters and topology; Helmholtz/reference calculations; phase and caloric primitives; exact explicit parameter partials; domain/applicability and artifact identity | reaction equilibrium, Ceres, application datasets, or promotion |
| Equilibrium | source-to-Provider reference transformation; reacting-liquid and reactive-bubble solves; exact implicit state sensitivities; solver/numerical/physical certificates | parameter selection, residual weights, regression objectives, or application promotion |
| Regression | source-bound observation validation; parameter transforms and sharing execution; residual/Jacobian assembly; Ceres; rank, conditioning, active-bound and confirmation diagnostics; immutable authority-neutral results | chemistry defaults, a second equilibrium formulation, copied EOS/reaction equations, or Provider catalog mutation |
| Validation campaign | installed-artifact black-box replay and durable cross-package evidence authored by the accountable package task | production algorithms or private source imports |

Once organization doctrine admits the proposed transport, the downstream
application constructs the model-bound evaluator from exact installed Provider
and Equilibrium artifacts and supplies the versioned process-local handle
described by Regression issue 15. Regression does not import or link
Equilibrium.

## 3. Current source contract

The authoritative MEA application baseline is
`MEA-Thermodynamics@c3a92740dfd1a53b5ef5197d26cc0af4d227afe3`.
Later exploratory SciPy fits are rejected evidence and do not define this
contract.

### 3.1 Chemistry

Canonical species order:

```text
CO2, MEA, H2O, MEAH+, MEACOO-, HCO3-, CO3^2-, H3O+, OH-
```

Canonical charges:

```text
0, 0, 0, +1, -1, -1, -2, +1, -1
```

Canonical reactions:

```text
2 H2O                <=> H3O+ + OH-
CO2 + 2 H2O          <=> HCO3- + H3O+
HCO3- + H2O          <=> CO3^2- + H3O+
MEACOO- + H2O        <=> HCO3- + MEA
MEAH+ + H2O          <=> MEA + H3O+
```

The application source contract is
`data/reference/MEA/manifests/chemical_reaction_source_contract.json`,
SHA-256
`77c6c0a705aa2f5bcc5f9ff90722ffadfe9b2ef0ae3cf41b8a1ac8fecbd286d8`.
Its Provider-basis transform is not yet complete. The source spelling
`CO3--` and runtime spelling `CO3^2-` must be reconciled by one explicit
identity mapping in the installed application contract, not by a Regression
alias.

### 3.2 Pressure and speciation partition

The frozen grouped split is
`data/reference/MEA/manifests/grouped_split_manifest.csv`, SHA-256
`af205ad5968667cf25dc9205d780738035769664a94cc9a421cd3c67148ff804`.

| Observation family | Training | Reserved |
|---|---:|---:|
| CO2 partial-pressure rows | 89 | 167 |
| Speciation states | 58 | 53 |
| Total states | 147 | 220 |
| Eligible speciation scalar records | 208 | 100 on 16 states |
| Maximum eligible pressure plus speciation residuals | 297 | 267 |

The reserved total is 267, not 435. The stale value counted 168 speciation
records whose manifest explicitly marks them ineligible.

Training memberships:

- pressure: Aronu `vle_obs_0034` through `vle_obs_0069`; Hilliard
  `vle_obs_0119` through `vle_obs_0142`; Idris `vle_obs_0162` through
  `vle_obs_0171`; Ma'mun `vle_obs_0246` through `vle_obs_0264`;
- speciation: Böttinger states 030 through 068; Matin states 001 through 019.

Reserved memberships:

- pressure: Aronu 0001–0033 and 0070–0106; Hilliard 0107–0118 and
  0143–0161; Jou 0186–0233; Xu 0265–0272 and 0279–0288;
- speciation: Böttinger 001–029 and Jakobsen 001–024. Only Jakobsen
  009–024 currently carries eligible quantitative validation observations.

The speciation membership file is
`data/reference/MEA/manifests/speciation_target_membership.csv`, SHA-256
`1e267adc62a6572d029caf92cfa39529cd38c9c0475fdf2f66a9546dc92c8f9b`.

### 3.3 Volumetric partition

The separately frozen volumetric campaign contains:

- 153 future training densities: 63 Amundsen, 56 Augusto, and 34 Dhage;
- 78 reserved densities: 40 Amundsen, 28 Augusto, and 10 Dhage.

The 44 derived excess-molar-volume values are diagnostics derived from parent
density rows. They must not be added as simultaneous independent residuals.
The volumetric campaign is a parameter-identification stage, not evidence that
the coupled pressure/speciation path is ready.

## 4. Observation contracts

### 4.1 State inputs

Every state must bind:

```text
row_id
source_id and immutable source hash
partition and group identity
temperature_K
pressure specification and pressure role
unloaded_MEA_mass_fraction
CO2_loading_mol_per_mol_MEA
apparent carbon, MEA, and water totals
ordered species/reaction/reference identities
installed artifact fingerprints
```

Loading, temperature, and unloaded MEA fraction are state specifications. They
are not residuals unless a future source explicitly reports them as uncertain
observations and admits a data-reconciliation contract.

### 4.2 CO2 partial pressure

Each pressure row must add:

```text
observed_p_CO2
unit
measurement_origin:
  direct_partial_pressure |
  calibration_derived_partial_pressure |
  total_pressure_derived |
  model_derived
measured_primitive
derivation and solvent-vapor correction
gas-analysis method and reported vapor composition
uncertainty and covariance when reported
source table/figure/row locator
residual transform, scale, and rationale
```

`model_derived` rows are never fit targets. The current source packet does not
yet classify pressure origin. Only `vle_obs_0171` presently carries a numeric
pressure uncertainty, 0.077 kPa; that makes it a useful tracer candidate but
does not freeze it.

### 4.3 Liquid speciation

Admitted roles are:

- positive species equality;
- source-defined linear aggregate equality, principally `MEA + MEAH+`; and
- one-sided upper censoring for a reported zero or below-detection result.

Training contains 136 positive equalities, 58 aggregates, and 14 censored
records. Reserved evidence contains 83 positive equalities, 16 aggregates,
and one censored record.

Aggregate coefficients, measurement identity, and covariance must be explicit
to avoid double counting. A reported zero without a numerical detection limit
cannot become `log(0)`, an invented zero equality, or an invented censor bound.
The current 14 training and one reserved zero records are therefore
non-executable until MEA acquires source-backed bounds or demotes them and
recomputes the executable observation counts and hashes. The 297/267 counts
remain maxima until that censor and pressure-origin adjudication finishes; the
grouped split hash changes only if state membership changes.
Balance-inferred, ambiguous, contextual, calibration-derived, and
`target_eligible=no` records do not enter Ceres.

### 4.4 Heat of absorption

Heat is `NOT_ADMITTED_NO_SOURCE_CONTRACT`. There are no qualified rows.
Before heat can join the objective, MEA must freeze:

- row identities and primary-source locators;
- differential, integral, or semi-differential definition;
- calorimetric versus derived measurement role;
- sign convention;
- initial and final loading/reference states;
- constant-pressure or other process path;
- energy unit and per-mole-absorbed-CO2 basis;
- uncertainty, covariance, partition, residual scale, and promotion gate.

Provider must then expose the required caloric primitives and exact temperature
and active-parameter partials, and Equilibrium must expose the corresponding
total sensitivity through the reacting state. This capability is outside
Equilibrium issues 36–38.

## 5. Thermodynamic and derivative contract

For a fixed application state and physical fitted parameters `theta`,
Equilibrium owns the reacting-liquid optimization and its complete active-set
KKT system:

```text
H(w, theta) =
  [gradient_z(L),
   c(z, theta),
   g_active(z, theta)] = 0
```

where `w` contains the reduced amount chart, `log(V)`, multipliers, and the
active-set variables required by the Equilibrium formulation; `c` is the
implemented equality-constraint block and `g_active` is the fixed active-bound
block. The equality block contains the compiled independent material balances;
the ionic amount chart enforces electroneutrality identically and charge is
independently certified after the solve. Reaction and pressure equilibrium
appear once through the implemented stationarity formulation—the fixed-pressure
objective contains the `P*V/(R*T)` term—and are not appended again as duplicate
constraints. On a fixed smooth active set:

```text
H_w * dw/dtheta = -H_theta
dw/dtheta = -H_w^{-1} H_theta
```

Equilibrium then extracts the state-coordinate rows `dz/dtheta`. For any
primitive observation `y(z, theta)`, the exact total derivative consumed by
Regression is:

```text
dy/dtheta = partial(y)/partial(theta)
          + partial(y)/partial(z) * dz/dtheta
```

Provider-parameter sensitivities require exact Provider mixed partials of the
state gradient with respect to every active physical parameter, exact
parameter derivatives of packing/domain functions, and exact derivatives of
the transformed reference vector. Missing terms make the Jacobian unavailable;
neither Equilibrium nor Regression substitutes finite differences or frozen
speciation.

CO2 pressure additionally requires Equilibrium to re-solve the reacting liquid
at each trial pressure, solve and certify the installed-Provider incipient
vapor condition, and root-find the reactive bubble boundary. A
frozen-speciation flash is initialization only. The complete derivative

```text
d(y_CO2 * P_bubble)/dtheta
```

also requires exact Provider parameter partials and implicit sensitivities for
the vapor-incidence residual and vapor-composition solve, together with the
source-reference pressure and parameter derivatives. A missing liquid, vapor,
reference, or incidence term makes the entire pressure Jacobian unavailable.
The bubble derivative also fails closed for vanishing or ill-conditioned
`partial(B)/partial(P)`, a singular vapor-composition Jacobian, root
swap/coalescence, or a vapor active-set/topology change. The liquid KKT
condition threshold does not certify this separate bubble-root system.

## 6. Equilibrium value and certificate contract

The future installed composed evaluator receipt must retain:

- `T` in K and `P` in Pa;
- ordered species IDs, charges, amounts in mol, mole fractions, and volume in
  m3;
- Provider component order, parameter fingerprint, reference identity, and
  artifact identity;
- ordered reaction and balance rows;
- source standard-state, activity-scale, conversion, and reaction-constant
  provenance;
- amount-chart topology and all active lower, upper, trace, and domain bounds;
- independent solver, numerical, physical, local-minimum, search, root,
  predictive, and globality statuses.

The current public homogeneous-liquid result already exposes the core state and
certificate values but does not expose all ordered balance/reaction matrices,
artifact hashes, sensitivity active sets/chart topology, or
search/root-completeness fields above. For homogeneous liquid results, search
and root fields are `not_applicable`; for reactive-bubble results they are
mandatory.

Current Equilibrium value gates independently recompute:

- balance and charge closure at `1e-9` absolute;
- pressure closure at `1e-8` relative;
- reaction affinity at `1e-7`;
- KKT stationarity and complementarity at `1e-7`;
- positivity, packing, Provider domain, and reduced-Hessian local-minimum
  status.

Sensitivity availability is separate from value validity. A Jacobian fails
closed for an uncertified primal state, singular system, infinity-norm
condition estimate above `1e6`, bound activity within the `1e-7` margin,
trace/structural-face activity, topology change, missing Provider mixed
partial, or missing transformed-reference derivative.

Equilibrium issue 36 is closed for exact conditioned sensitivities to balance
totals, Provider-basis `ln(K)`, and pressure. It does not yet supply active
Provider-parameter sensitivities or source-bound reference derivatives.
Issue 37 remains open for the reactive bubble. Issue 38 remains open for
installed-artifact MEA liquid and bubble evidence.

Current public `chemical_equilibrium` execution exposes values only.
Sensitivity data remain in a private underscored payload and source-bound
sensitivities fail closed. Before Regression issue 15 can bind the composed
evaluator, Equilibrium must expose one typed, versioned installed
value/Jacobian contract—not the private `_chemical_equilibrium` seam—that binds
wheel and ABI identity, state and parameter ordering, units, reference
identity, chart topology, and Provider fingerprint.

## 7. Regression residual and Jacobian assembly

Regression reuses the existing ordered `[N fitted | Q lifted]` Ceres core,
result family, native module, and native target. Equilibrium eliminates the
reacting state, so the coupled MEA path has `Q = 0`.

The application chooses one admitted transform per observation:

```text
positive equality: identity or log transform
linear aggregate:  a^T y
upper/lower censor: frozen differentiable one-sided rule
```

The log base, scale, family normalization, and any uncertainty interpretation
remain MEA inputs. Equal weights of 1.0 and the historical regularization
weight 0.003 are not source-backed defaults and shall not be inferred.

For the evaluator's ordered primitive values `y`, exact physical-parameter
Jacobian `Y_theta`, sharing map `S`, and affine solver transform `D`,
Regression assembles:

```text
J = W * h_y * Y_theta * S * D
```

Every requested fitted column must be exact, finite, ordered, and complete.
Censored rows that are locally inactive contribute no local rank.

## 8. Parameter-identification stages

The old 12-parameter pressure/speciation block is not the production plan.

| Stage | Active coordinates | Current role |
|---|---|---|
| Shared-cation analog | `MEAH+::sigma`, `MEAH+::epsilon_over_k` | identify a common cation contribution from admitted direct analog densities after immutable counterion parameters/covariance or declared joint regularization exist |
| Carbamate correction | `MEACOO-::sigma` | smallest correction from electroneutral reactive-liquid evidence after an immutable transfer prior and scale exist; isolated-ion states and directional Maiti context are not residuals |
| Reactive refinement | the three coordinates above | frozen objective is Amundsen density plus eligible speciation; pressure is a rejection guard unless a pre-fit amendment freezes metrology, scaling, sensitivities, and the revised objective |

Initial safeguard contract:

| Parameter | Historical start | Bounds |
|---|---:|---:|
| `MEAH+::sigma` | 3.48508556586 angstrom | 2.0–5.8 angstrom |
| `MEAH+::epsilon_over_k` | 232.687201645 K | 50–950 K |
| `MEACOO-::sigma` | 3.53543525721 angstrom | 2.0–5.8 angstrom |

These starts are provisional seeds, not scientific authority. Affine coordinate
scales remain unfrozen. Segment counts, Born diameters,
`MEACOO-::epsilon_over_k`, reaction constants, trace-ion parameters, and all
binary interactions remain fixed initially.

`k_ij(MEAH+,MEACOO-)` may replace, not augment, the carbamate correction only
after a preregistered sensitivity amendment. A quadratic transfer prior for
`MEACOO-::sigma` remains blocked on an immutable uncertainty scale and generic
regularization support.

Frozen density scaling is application input:

- Augusto uses the greater of the row standard deviation and reported relative
  uncertainty times density;
- Dhage uses 0.0005 g/cm3;
- Amundsen uses 0.0005 g/cm3 for unloaded rows and 0.002 g/cm3 for loaded rows;
- each source/salt family is normalized by its admitted row count.

## 9. Executable sequence

### Gate 0 — application contract

MEA must:

- complete the Provider-basis reaction transformation and carbonate identity;
- enforce `target_eligible`;
- propagate eligible aggregates into the executable payload;
- acquire source-backed detection bounds for censored rows or demote them, then
  recompute executable counts and observation hashes;
- classify pressure metrology;
- freeze source-backed row pressure or an explicitly bounded pressure
  convention and a 2-oxazolidone inclusion/exclusion rule before selecting the
  sentinel or tracer state;
- freeze affine and residual scales, starts, bounds, sharing, and exact tracer
  rows;
- remove heat from the active tracer or admit a complete calorimetric packet.

### Gate 1 — installed reacting-liquid value sentinel

Equilibrium issue 38 first retains one source-complete fixed-`T,P` aqueous MEA
liquid state using the nine-species/five-reaction bundle and exact installed
Provider artifact. It must return eligible speciation primitives and every
liquid certificate. This is value evidence, not a fit.

### Gate 2 — exact active-parameter sensitivities

Provider supplies the exact explicit partials for only the preregistered
coordinate subset. Equilibrium extends issue 36's implicit solve to those
columns and the source-reference transform. Value-valid but
Jacobian-unavailable states remain unusable for Ceres.

### Gate 3 — reactive pressure

Equilibrium issue 37 supplies the certified reactive-bubble result and exact
total derivatives needed for `p_CO2 = y_CO2 P_bubble`. Issue 38 then retains
installed MEA bubble evidence.

### Gate 4 — cross-package evaluator transport

Equilibrium first exposes the typed installed public value/Jacobian contract in
section 6. After the organization-level transport decision, Regression issue
15 binds a downstream-composed, model-bound, versioned evaluator through the
admitted process-local transport. It supports value-only and
value-plus-Jacobian responses, deterministic ordering, batching, complete
fingerprints, per-row failures, and no Regression-to-Equilibrium dependency.

### Gate 5 — reduced coupled tracer

Regression issue 16 is corrected to fit one or two predeclared coordinates
against exactly:

- one admitted CO2 partial-pressure equality; and
- one eligible speciation equality or aggregate.

The two-row sensitivity matrix must have rank `N` at every declared start
before optimization. If no source/physics-selected subset with `N <= 2` has
that rank, the tracer remains blocked; rows, weights, priors, or parameters are
not changed after observing the result.

### Gate 6 — staged identification and reactive refinement

Run the shared-cation analog and carbamate-correction stages only through exact
installed derivatives and their frozen source contracts. Then fit the
three-coordinate reactive block against the preregistered training
observations. Pressure may enter the objective only after metrology and scaling
are frozen; otherwise it remains a coupled rejection guard.

The full staged campaign uses 32 deterministic scrambled-Sobol starts with seed
`390035`, plus the provisional historical seed. Materially distinct best
basins, incomplete rank without profile-supported identifiability, material
start/prior/bound dependence, or a safeguard-bound solution rejects promotion.

### Gate 7 — reserved evaluation

Evaluate without refitting:

- the 220 reserved pressure/speciation states, reporting 37 context-only
  speciation states separately from the 16 quantitatively eligible states; and
- the 78 reserved volumetric rows.

Do not combine derived excess volumes with parent density residuals. Preserve
every failed, skipped, inapplicable, and evaluated row.

Any failed eligible training state rejects the fit. Every failed reserved state
counts as a failed prediction and cannot be relabeled skipped or inapplicable
after execution.

### Gate 8 — later heat expansion

Only after the source, caloric Provider, Equilibrium sensitivity, and
observation-scale gates in section 4.4 are complete may heat join the same
parameter vector and Ceres objective. Its addition requires a new rank and
conditioning preflight; it does not reinterpret the earlier tracer.

## 10. Performance and execution shape

One evaluator batch groups observations only when their complete
specifications and equilibrium operations are identical:

1. perform one top-level evaluation per logical state and operation;
2. reuse its certified values for every attached pressure, speciation,
   aggregate, censor, density, or later heat observation;
3. request a residual-only response, without returned fitted-parameter
   columns, for Ceres residual-only calls;
4. request all exact fitted columns once for Jacobian calls;
5. return deterministic row order and one status for every input.

A fixed-`T,P` liquid, reactive-bubble boundary, and calorimetric path are
different operations even when temperature, loading, and composition match.
Every bubble pressure trial still re-solves the reacting liquid. A
residual-only response may omit returned fitted-parameter columns, but
Equilibrium still consumes the exact state derivatives and Provider tensors
required internally by Ipopt and certification; it cannot replace them with a
value-only phase endpoint.

No nested solver parallelism is allowed. Bounded concurrency requires an
explicit installed reentrancy capability. Per-fit budgets and cancellation
must be task-owned. The historical IDAES/eNRTL workflow is architectural
corroboration only; it is not copied, imported, or used as a runtime
dependency.

## 11. Result and acceptance semantics

Every result separately reports:

- upstream Equilibrium solver, numerical, physical, local-minimum, search,
  root-completeness, and globality status;
- Ceres termination and Regression numerical validity;
- full Jacobian rank, conditioning, active bounds, KKT/bound diagnostics, and
  primary/confirmation-start agreement;
- pressure, each speciation species/aggregate/censor, density, and later heat
  residual summaries separately;
- exact input, evaluated, skipped, failed, training, and reserved accounting;
- workflow validity, scientific comparison, and predictive status.

Optimization success never converts a rejected or uncertified state into a
valid row. A density improvement accompanied by degraded reserved pressure or
speciation fails. Current legacy median-error gates are application policies,
not experimental-uncertainty limits, and require explicit preregistration
before reuse.

A promotion candidate additionally retains one complete successor receipt:
exact Provider, Equilibrium, and Regression versions, commits/trees, installed
wheel or ABI hashes, and capability fingerprints; ordered parameter identities,
units, values, bounds, sources, transforms, and fitted domain; observation and
split hashes; reaction/reference identities; all multistart diagnostics;
Equilibrium certificates; complete training/reserved residual records; and
final bundle fingerprint.

No result establishes global parameter uniqueness, uncertainty, prediction,
Provider-catalog authority, or promotion by itself.

## 12. Current readiness

| Dependency | Status | Smallest missing work |
|---|---|---|
| Regression multi-parameter Ceres core | ready for admitted exact evaluators | none |
| Generic observation equations | design frozen | bind real downstream evaluator in issue 15 |
| MEA chemistry/source transform | blocked | complete Provider-basis transform and explicit carbonate identity |
| MEA pressure metrology | blocked | classify every candidate pressure row and exclude model-derived rows |
| MEA speciation payload | blocked | enforce eligibility and executable aggregate/censor semantics |
| MEA heat | not admitted | source-complete packet plus caloric value/sensitivity contract |
| Equilibrium reacting-liquid values | implemented generically, not yet MEA-evidenced | one installed source-complete MEA sentinel |
| Equilibrium active-Provider-parameter sensitivities | blocked | exact Provider mixed partials and source-reference derivatives |
| Equilibrium reactive bubble | blocked | issue 37 |
| Installed MEA campaign | blocked | issue 38 after its liquid and bubble gates |
| Regression downstream transport | blocked | issue 15 after exact installed evaluator artifact |
| Reduced coupled fit | blocked | corrected issue 16 after gates 0–4 |
| Full mixed campaign | not ready | successful tracer, frozen three-coordinate contract, and rank-sufficient admitted rows |

The next Regression-side work that does not depend on Equilibrium is the
generic strict dataset/control and preparation/preflight work in issues 21 and
20. It must remain chemistry-free and must not fabricate missing MEA
observations or derivatives.

## 13. Required issue reconciliation

- Regression issue 16 must replace “CO2-loading plus heat” with “CO2 partial
  pressure plus eligible speciation,” correct the reserved maximum from 435 to
  267, and retain `N <= 2`.
- Regression issue 15 remains the transport gate; it must not be bypassed by a
  Python callback or source checkout.
- MEA issues 12 and 13 must bind the corrected split, observation roles,
  three-coordinate plan, and exact installed artifact identities.
- MEA issues 35–37 own the remaining source, Provider-basis, and sentinel
  preparation.
- Equilibrium issue 36 is complete only for its declared totals, `ln(K)`, and
  pressure sensitivity scope. Issues 37 and 38 remain open.

This document controls the coordinated sequence. Application manifests,
installed artifact receipts, and owner-specific implementation contracts
remain authoritative for their own content.
