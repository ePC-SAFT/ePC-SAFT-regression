# Figiel Current-Catalog Born-Diameter Tracer Design

Status: authority-neutral design ready for permanent-lab review; runtime,
Provider implementation, fitting, artifacts, and authority remain unauthorized.

Migration decision D-023 authorizes this documentation-only checkpoint at
Migration commit `e0d6930228ad8d9f7c98c419110bbb01bde3a62f`, tree
`582c4c02bb864069e379c9fd2373cf8887ffc1c6`. This document is the sole
Regression design owner for the first Figiel tracer. It creates no second plan
or generic parameter-family contract.

## Decision and problem shape

The first tracer is one bounded Ceres problem with one contiguous five-vector
parameter block, five dimensionless variables, and five dimensionless
residuals in fixed order:

```text
(Li+, Na+, K+, Cl-, Br-).
```

Each variable activates only its ion's Born diameter in one Provider-bound
model. Each residual compares the corresponding single-ion water solvation
Gibbs energy with one SI Table S5 reported average. The problem is separable,
but one `5 x 5` solve is the minimum package shape because it produces one
result, one solver execution, and one explicit system-rank adjudication.

Rejected alternatives are:

- five independent scalar solves, which multiply orchestration and result
  semantics without adding evidence;
- an analytic or bracketed root solver, which creates a second fitting route
  instead of using the repository's admitted Ceres owner; and
- a generic active-parameter registry or overlay, which admits families and
  identity policy not required by these five targets.

The eventual implementation may add one family-specific cost owner and one
family-specific immutable result to the existing `_native` module and native
target. It must not generalize the accepted pure-saturation cost or result into
a universal regression framework. Methane and ethane behavior remains
unchanged.

## Immutable source contract

The approved source subject is Validation commit
`8944d34f7002cda1bb8760e606cc1f11696f58cd`, tree
`6c8fd350dcd6bfdd7be1918f73fd33a23e2070dd`:

| Record | SHA-256 |
|---|---|
| `data/figiel-2025-regression-target-ledger.csv` | `f405a3e48d21cd979a8dd480d5f8cb3be40754f5d6babf368b505b5f305607f0` |
| `data/figiel-2025-parameter-provenance.csv` | `932e8baa90fcefbaa8c3a8730cdeadd83a4c01f0a3b109f4e4cd0319aee9312b` |
| `data/figiel-2025-regression-target-ledger.yaml` | `8ea06c6ca5452d01448a03f9a76cf7d0c35bb99c9abe23ccb1729d56c71d468f` |

The literature identity is Figiel, Yu, and Held (2025), main-paper DOI
`10.1021/acs.iecr.5c00475` and official supporting-information DOI
`10.1021/acs.iecr.5c00475.s001`. The controlling numerical locator is SI Table
S5, PDF page 9 of 10, `lit` reported-average column. The official SI PDF has
SHA-256 `005b38ed566ec3c09b87e1ca3a9dd6eeafc9ba75e1a30b9322291d770bb93895`.
The byte-identical retained main-paper Markdown has SHA-256
`ce80533925a91bc59d8d0d8056113c40611ca26c2edf04aced76986d50bd4bae`;
local paths are read-only execution metadata, not package authority.

Exactly these five ledger rows are residual targets:

| Target ID | Ion | Reported target [kJ mol^-1] | Exact runtime value [J mol^-1] | Published Table 3 Born diameter [angstrom] |
|---|---|---:|---:|---:|
| `figiel2025-s5-Lip-reported-average` | `Li+` | -486.2 | -486200 | 2.784 |
| `figiel2025-s5-Nap-reported-average` | `Na+` | -381.1 | -381100 | 3.445 |
| `figiel2025-s5-Kp-reported-average` | `K+` | -309.1 | -309100 | 4.150 |
| `figiel2025-s5-Clm-reported-average` | `Cl-` | -314.9 | -314900 | 4.100 |
| `figiel2025-s5-Brm-reported-average` | `Br-` | -290.9 | -290900 | 4.480 |

Conversion from kJ/mol to J/mol is exact multiplication by 1000. The sign is
preserved: negative is favorable gas-to-solution transfer. The source basis is
the x-process at `298.15 K` and `1 bar`, from gas at 1 bar to a hypothetical
dilute-ideal solution at 1 bar with `x_ion -> 0`. A reported average is used as
reported and is not recomputed from its literature values.

The other 27 current-ion SI Table S5 rows are provenance/support only. All 164
Hamer--Wu MIAC rows belong to later `f_k` and `k_ij` stages. No target has a
pointwise uncertainty or approved comparison threshold, and no such value is
inferred here.

## Approved Provider boundary

The Provider prerequisite is the permanent-lab-approved design at commit
`da9660481f08bb5557cc03da528edec15cc15e1f`, tree
`e34575ae646c40f3fb63b5994c957e30bb035f69`. It specifies one future
model-bound callback, `evaluate_ion_solvation_born`, that returns:

```text
G_i(d_i)  = Delta_solv,x G_i^infinity             [J mol^-1]
G'_i(d_i) = d G_i / d d_i at fixed T, P, path     [J mol^-1 angstrom^-1].
```

The callback fixes `T = 298.15 K`, `P = 100000 Pa`, water, the four-molality
reference sequence `(1e-6, 1e-8, 1e-10, 1e-12) mol kg^-1`, and all catalog
records except the active ion's trial Born diameter. It returns the terminal
`1e-12 mol kg^-1` value only when the maximum final change in either ion's
`ln(phi)` is at most `5e-5`. The Provider owns the EOS, association state,
density roots, reference composition sequence, convergence adjudication,
terminal fugacity selection, and the exact total fixed-pressure derivative.
Regression copies none of those equations or algorithms.

The future callback is design-only; no installed header, wheel, or runtime
entry exists. Runtime work therefore remains dependency-gated after this
design review.

The admitted models and returned fingerprints are:

| Active ion | Bound component order | Required source fingerprint |
|---|---|---|
| `lithium-cation` | `(water, lithium-cation, chloride-anion)` | `sha256:1bb528ebe8f5612757e148608fc55821f9fb03737dbcec6d0bc4fffd0f4cbc4c` |
| `sodium-cation` | `(water, sodium-cation, chloride-anion)` | `sha256:7c637771bc9f717b8f47b44bb2a61044c3fe83084dca7c3c16102fba0989912d` |
| `potassium-cation` | `(water, potassium-cation, chloride-anion)` | `sha256:d29cef0c0f63034436d547d0aafa57934effe06783c8dffd89c94fa85e6940f4` |
| `chloride-anion` | `(water, chloride-anion, sodium-cation)` | `sha256:7551f1eee5903b66061cf7520f3bb59b169896ce372f3df3d48aa7ec778c39d4` |
| `bromide-anion` | `(water, bromide-anion, sodium-cation)` | `sha256:70ae04599dfa8338175e793bac6b9e4dfad37a9b96a568b5484dc87f104ef1a9` |

Water is component 0, the active ion is component 1, and its fixed counterion
is component 2. Every model must have exactly three components, the
`ionic-ssm-ds` topology, water as its only neutral solvent, and charge order
`(0,+1,-1)` or `(0,-1,+1)`.

The ordered five-model tuple is the workflow input. Regression rejects the
wrong count, order, component topology, fingerprint, callback status parity,
terminal molality, or reference convergence before interpreting a fit. A
missing callback or incompatible table/result size fails before Ceres through
the existing runtime-error path with reason `provider capability unavailable`;
there is no compatibility fallback or new error-code vocabulary.

## Variables, bounds, scales, and starts

For each ion, the dimensionless Ceres coordinate is

```text
d_i = 3.0 angstrom + (1.0 angstrom) z_i.
```

Every diameter has the closed engineering bound

```text
1.0 angstrom <= d_i <= 6.0 angstrom,
-2 <= z_i <= 3.
```

This is not a universal ion-size claim. It is a common whole-angstrom
engineering enclosure of all ten published Table 3 Born diameters, whose
verified range is `1.218` to `4.985 angstrom`. The lower bound is the largest
positive whole angstrom below the source minimum; the upper bound leaves more
than one angstrom above the source maximum. Bounds must not be widened after
seeing a result.

Published fitted diameters are comparison anchors only. They are never starts
or hidden initial values. The deterministic start schedule is:

| Solve | Diameter vector [angstrom] |
|---|---|
| primary | `(3, 3, 3, 3, 3)` |
| lower confirmation | `(2, 2, 2, 2, 2)` |
| upper confirmation | `(5, 5, 5, 5, 5)` |

The common starts are source-independent length values inside the predeclared
bound. No random seed, continuation, or per-ion start tuning is allowed.

## Residuals and exact Jacobian

For target `G_i^target` in J/mol, define the magnitude scale

```text
s_i = abs(G_i^target) [J mol^-1]
r_i(z_i) = (G_i(d_i) - G_i^target) / s_i.
```

All five residuals are dimensionless. Target-magnitude scaling keeps the five
reported-average equations comparable; it is not an uncertainty weight. No
robust loss is used. Ceres minimizes

```text
1/2 sum_i r_i^2.
```

The exact Jacobian consumed from the approved first-derivative callback is

```text
J_ij = d r_i / d z_j
     = delta_ij * G'_i(d_i) * (1 angstrom) / s_i.
```

Thus `J` is exactly diagonal in the frozen ion order. The zero off-diagonal
entries follow from the five independently bound Provider models, not from a
copied thermodynamic equation. No second derivative, numerical production
derivative, density-root derivative, or equilibrium dependency is consumed.

## Ceres and local sensitivity contract

The future implementation reuses the accepted deterministic Ceres policy:

```text
linear solver       = DENSE_QR
threads             = 1
logging             = SILENT
max iterations      = 500
function tolerance  = 1e-10
gradient tolerance  = 1e-10
parameter tolerance = 1e-10
```

Each evaluation retains raw `G_i`, raw `G'_i`, raw error in J/mol, scaled
residual, scaled diagonal Jacobian entry, callback diagnostics, and fingerprint.
The final `5 x 5` Jacobian is decomposed by SVD. Reusing the existing package
rank rule, singular value `sigma_k` is accepted when

```text
sigma_k > sigma_max * max(5, 5) * epsilon_double * 100.
```

Full parameter rank must be exactly 5 for the primary and both confirmation
solutions. The result retains every singular value, the threshold,
`sigma_min`, `sigma_max`, and `kappa_2 = sigma_max/sigma_min`. No separate
condition-number cutoff is invented: nonfinite conditioning or loss of rank
fails, while a finite full-rank condition number remains descriptive. Because
the matrix is diagonal, the result also identifies the least-sensitive ion.
Rank 5 supports only local linear identifiability of these five independent
coordinates. It does not establish global uniqueness or global
identifiability.

For scaled bounds `z_L=-2` and `z_U=3`, retain each final bound distance. A
parameter is numerically active when

```text
min(z_i - z_L, z_U - z_i)
  <= sqrt(epsilon_double) * max(1, abs(z_L), abs(z_U)).
```

An active bound rejects the recovery claim. This is a binary64 resolution
test, not a physical confidence interval.

## Confirmation and status semantics

All three declared starts must run without changing bounds, scales, targets,
or solver controls. Primary, lower, and upper solutions must each have Ceres
`CONVERGENCE`, a usable finite solution, finite complete Jacobian columns,
strict cost reduction, full rank 5, parameters inside and inactive at bounds,
and callback success for all five models.

Numerical confirmation requires, for each confirmation relative to primary,

```text
max_i abs((d_i^confirmation - d_i^primary) / 1 angstrom) <= 1e-5
max_i abs(r_i) <= 1e-8 for every solve.
```

Initial/final costs and symmetric relative cost deltas are retained but are not
an acceptance gate because an exactly determined root problem can drive both
costs into the binary64 floor. A failed start is reported; it is not replaced.

The immutable result keeps four meanings separate:

1. `solver_converged`: the primary Ceres termination, usability, finiteness,
   rank, bounds, cost-reduction, and callback gates pass.
2. `numerically_converged`: solver convergence plus both confirmation and
   `1e-8` scaled-residual gates pass. The `1e-8` value is the repository's
   existing dimensionless numerical-closure policy, not source uncertainty.
3. `physically_valid`: numerical convergence plus exact source/model/state/
   x-process identities, Provider reference convergence, and both recovery
   comparisons below pass. This status means source-bound workflow recovery,
   not general physical or predictive validation.
4. `predictive_status` remains
   `NOT_ADJUDICATED_NO_APPROVED_HELD_OUT_CUTOFF`. All five targets are training
   equations; no held-out observable or predictive cutoff exists.

The source reports targets to `0.1 kJ/mol` and diameters to
`0.001 angstrom`. Therefore the two recovery comparisons are fixed before any
fit as decimal-reporting round-trip gates:

```text
abs(G_i^fit - G_i^target) <= 50 J mol^-1
abs(d_i^fit - d_i^published) <= 0.0005 angstrom.
```

These half-last-digit values express source resolution only. They are neither
experimental uncertainty nor statistical confidence intervals. No tolerance
may be changed after observing a fit; failure returns a failed recovery claim.

## Singular ownership and future result

This design changes no runtime or public interface. If separately authorized,
the minimum family surface is exactly:

- one immutable closed specification, `FIGIEL_BORN_DIAMETER_TRACER_V1`;
- one workflow, `fit_figiel_born_diameters`, accepting the exact ordered
  five-model tuple; and
- one immutable `BornDiameterFitResult` containing all five parameters,
  residuals, sensitivities, ranks, confirmations, statuses, and reasons.

The existing `records.py` remains the source/specification owner and
`workflow.py` remains the result/status owner. The single family-specific
native cost owner is `src/epcsaft_regression/native/born_diameter_fit.cpp`
inside the existing `epcsaft_regression._native` module and its single CMake
target; it must not be added to the already 920-line pure-saturation cost file.
No per-ion class, target registry, backend selector, second native module, or
second native target is admitted.

The future result binds the three Validation hashes, Provider commit/header/
wheel identities, five expected and observed fingerprints, target IDs and
values, state and convention, starts/bounds/scales, Ceres controls, costs and
iterations, raw and scaled residuals, raw derivatives, SVD/rank/conditioning,
active-bound distances, published-parameter deltas, confirmations, four
statuses, and ordered failure reasons. It contains no covariance, uncertainty,
catalog-write request, or global-identifiability field.

## Installed-artifact evidence plan

Runtime implementation cannot begin from this design alone. After permanent-
lab approval, Provider must first implement its approved single-owner callback,
return a stable commit plus exact wheel/header hashes, and pass its design's
value, derivative, reference-owner parity, ABI, failure, and negative-space
evidence. Migration must bind that exact artifact before Regression consumes
it.

A later Regression implementation must then:

1. retain only the compact five-row source record plus durable DOI/table
   locators and the three approved Validation hashes; it must not copy the
   407-row ledger or depend on Validation at runtime;
2. bind the installed Provider header and wheel only, reject a missing or
   truncated callback, and preserve one native module and target;
3. check the residual/Jacobian block against centered callback-value
   differences with step halving at non-published trial diameters, using the
   Provider design's approved five-point criterion before exact residual
   scaling:

   ```text
   abs(G'_CppAD - G'_FD)
     <= max(1e-3 J mol^-1 angstrom^-1,
            20 abs(G'_h - G'_(h/2)),
            2e-8 abs(G'_CppAD));
   ```

   numerical differences remain test-only;
4. run the three frozen starts and retain every gate and comparison above;
5. freeze methane and ethane accepted numerical parity and prove no new
   dependency, export alias, target registry, Provider source, sibling source,
   or copied Born/EOS/reference-sequence logic; and
6. build one immutable Regression wheel and prove wheel members, import origin,
   linkage, native-module/target counts, and absence of lab, Migration,
   Validation, or sibling-source runtime paths.

Validation later installs the exact Provider and Regression wheels in an
isolated environment. Its one bounded campaign independently binds the five
ledger targets and published diameters, checks sign/unit/x-process/state and
source-resolution round trips, replays all three starts, verifies fingerprints
and negative controls, and reports solver, numerical, physical/workflow, and
predictive statuses separately. It does not import package source or generate
production expected values. A passed campaign remains authority-neutral until
a later receipt and exact user approval.

Material negative controls include wrong target sign or kJ/J conversion,
changed target/order/hash, duplicated or missing ion, wrong component order or
fingerprint, unlisted H+/I-/sulfate/V models, neutral or organic solvent,
additional ions, wrong temperature or pressure, nonfinite/nonpositive trial
diameter, active bounds, callback ABI/status mismatch, reference-limit failure,
rank loss, and perturbed-start disagreement.

## Explicit negative space

This tracer does not admit `f_solv`, `k_ij`, dielectric suppression, neutral-
solvent parameters, ion `sigma` or `epsilon/k`, association parameters,
expanded ions or vanadium species, a generic overlay or registry, a combined
Tables 2–5 solve, Provider catalog persistence, a second Provider derivative,
another native target, a compatibility shim, prediction, uncertainty, global
uniqueness, global identifiability, or a runtime dependency on the lab,
Migration, Validation, or sibling source. Propane, binary-`k_ij`, and pressure-
resolution subjects remain immutable deferred provenance and are not gates.
