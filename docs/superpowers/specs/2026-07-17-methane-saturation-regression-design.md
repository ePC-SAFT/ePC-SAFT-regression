# Pure-Methane Saturation Regression Design

Status: stage-owned implementation design
Date: 2026-07-17
Authority: the bounded user approval delegated from task `019f5dba-f699-79c2-94ca-93fdf49d9b4d`

## Capability and ownership

This slice fits exactly three pure-methane PC-SAFT parameters—segment count
`m`, segment diameter `sigma [angstrom]`, and dispersion energy
`epsilon/k [K]`—to four retained NIST saturation rows. Regression owns the
immutable target records, parameter transformations, residual assembly, Ceres
solve, diagnostics, and candidate receipt. The installed `epcsaft` provider
wheel remains the only owner of Helmholtz evaluation and nonlinear
derivatives. The regression build consumes only its declaration-only public
header, and runtime calls only the model-bound `epcsaft.native_sdk.v1`
capsule.

This is a candidate capability for manager and validation review. It does not
transfer runtime authority, write parameters to a provider catalog, or claim
predictive admission.

## Source and target contract

The compact owned dataset is the nine-row NIST Chemistry WebBook methane
saturation table from 100 K through 180 K in 10 K increments. It retains the
exact NIST query, retrieval date, original units, and SHA-256 of the distilled
CSV. The four training rows are exactly 110, 130, 150, and 170 K; all nine rows
are reporting evidence.

Input construction rejects a changed header, missing unit declaration or
source identity, non-methane species, nonfinite or nonpositive values,
duplicate row IDs, duplicate temperatures, out-of-order temperatures, an
unexpected training partition, or an unsupported target family. No
composition is accepted or normalized: this route is pure methane only.

## Variables and transformations

The global physical parameter tuple is

```text
p = (m, sigma_angstrom, epsilon_over_k_kelvin)
start = (1.08, 3.555744, 157.5315)
lower = (0.5, 2.0, 50.0)
upper = (3.5, 5.0, 400.0)
scale = (0.1, 0.1, 10.0)
p_j = start_j + scale_j z_j
```

Ceres optimizes the dimensionless `z_j`; physical bounds are transformed
exactly to bounds on `z_j`. Each training row fixes `n = 1 mol` and owns two
dimensionless log-volume variables:

```text
V_liquid = V_liquid_start exp(u_liquid)
V_vapor  = V_vapor_start  exp(u_vapor)
V_liquid_start = M_methane / rho_liquid_observed
V_vapor_start  = R T / P_observed
```

`M_methane = 0.016043 kg/mol` and
`R = 8.31446261815324 J/(mol K)`. Static physical bounds
`V_liquid in [2e-5, 1e-4] m3/mol` and
`V_vapor in [1.5e-4, 0.1] m3/mol` keep both volumes positive and strictly
ordered. Their log bounds are computed from each row's start. The native
Jacobian applies `dV/du = V` and `dp/dz = scale` exactly.

The start is the retained M5 displaced methane tuple
`(1.0*1.08, 3.7039*0.96, 150.03*1.05)`. The molar mass is the retained
M5 formulation constant. Both are candidate inputs, not newly fitted data.

## Training residuals and exact Jacobian

For each training row, the provider is evaluated twice at the trial parameter
tuple, once at each row-local volume. The ordered raw residuals are

```text
P_liquid - P_observed                         [Pa]
P_vapor  - P_observed                         [Pa]
mu_liquid/(RT) - mu_vapor/(RT)                [dimensionless]
M_methane / V_liquid - rho_liquid_observed    [kg/m3]
```

Pressure residual scales are that row's observed pressure, the chemical-
potential scale is `1`, and the density scale is that row's observed liquid
density. Every residual has weight `0.25`, so the residual passed to Ceres is
`sqrt(0.25) * raw / scale`. These weights normalize the four residuals per row;
they are not uncertainty or accuracy claims.

With provider coordinate order `(n, V, m, sigma, epsilon/k)`, the exact
physical-coordinate derivatives are

```text
dP/dq_j       = -R T Phi_(V,q_j)
d(mu/RT)/dq_j = Phi_(n,q_j)
d(M/V)/dV     = -M/V^2.
```

The native cost function assembles all 16 by 11 Jacobian entries from the
provider Hessian and the exact transform chain rule. No finite difference or
alternate derivative path is present in runtime code. A test-only centered
directional difference checks `J v` against residual differences.

## Acceptance and diagnostics

Solver convergence is true only for Ceres `CONVERGENCE`, a usable finite
solution, finite initial and final cost, complete finite Jacobian columns, and
respected physical bounds. Numerical convergence is a separate confirmation
solve from perturbed row-volume starts; it requires both accepted solves and
agreement in scaled parameters and final cost at declared tolerances.
The confirmation starts multiply every liquid reference volume by `1.01` and
every vapor reference volume by `0.98`; agreement limits are `1e-5` in the
maximum scaled parameter coordinate and `1e-8` in relative final cost. Both
solves use Ceres `DENSE_QR`, one thread, silent logging, at most 500 iterations,
and function/gradient/parameter tolerances of `1e-10`.

Physical validity is separate again. Every final phase state must have finite
provider output, `Phi_VV > 0` (equivalently positive `dP/d(rho_m)`), liquid
volume below vapor volume, and relative phase-volume separation above `1e-3`.
Any provider error, unstable phase, coalescence, topology loss, or nonfinite
diagnostic produces an explicit failure reason.

The result retains start/final parameters, movement, bounds, active bounds,
full-Jacobian and parameter-column singular values/rank/condition diagnostics,
per-row raw and scaled residuals, phase volumes/densities/stability slopes,
Ceres termination and iteration counts, the confirmation-solve comparison,
and all failure reasons.

## Reporting evidence

At the final parameters, a fixed-parameter native Ceres solve predicts
saturation pressure and both phase volumes independently for each 100–180 K
reporting row. It uses the same pressure and chemical-potential relations and
exact provider-Hessian Jacobian. Observed values seed the solve but are not
residual targets. Reporting records retain predicted pressure, predicted
liquid density, both phase molar and mass densities, raw equilibrium residuals,
phase separation, stability, and
whether the row was in training.

Reporting uses `P_report = P_observed*exp(u_pressure)` with physical pressure
bounds `[1e3, 1e7] Pa`. A reporting row is physically valid only when each
pressure closure divided by observed pressure is at most `1e-8` in magnitude
and `|mu_liquid/(RT)-mu_vapor/(RT)| <= 1e-8`, in addition to solver, topology,
stability, and finiteness gates. These are numerical closure thresholds, not
uncertainty claims.

Held-out evidence is descriptive candidate evidence only. No accuracy cutoff
is invented; the receipt reports errors and keeps predictive admission owned
by validation.

## Package and negative space

The distribution contains one Python package, one CPython extension, and one
CMake native target linked to system Ceres. It has one public workflow
function and immutable typed input/result records. There is no binary `k_ij`,
association, electrolyte, reactive, generic registry, backend selector,
persistence layer, provider source, copied equation, provider private import,
or sibling-source build path.

The capsule consumer checks the capsule name and ABI prefix first. It reads the
tail only after `table_size` covers the parameterized fields, then checks
`parameterized_result_size` and the non-null function pointer before any
provider call. The extension retains the capsule reference for the complete
native solve.

Before Ceres construction, the native contract parser retains and exactly
validates the dataset/source IDs, source URL/citation/locator and hashes, row
IDs and values, units, parameter/residual names, transforms, bounds, fixed
amount, provider SDK ABI/table contract and source fingerprint, and solver
controls. Native row IDs and the complete compiled identity are returned and
checked by Python; results are never assigned source labels after the solve.
The public runtime cannot prove an artifact commit or wheel hash and therefore
does not claim either. The isolated candidate runner hashes and byte-compares
the installed artifacts and records that evidence in the canonical receipt.

## Verification and evidence boundary

Three durable tests cover the strict dataset/API contract, the installed-wheel
transport and native fit, and exact-Jacobian plus real-data anchors. Source and
isolated-wheel verification build against a provider installed from the pinned
wheel. Static scans reject sibling paths, private provider headers/imports,
provider sources, copied EOS terms, duplicate provider symbols, and excluded
target families. A candidate capability/architecture record and fit receipt
are generated for review but cannot accept themselves. The existing private
native test surface exposes the training residual/Jacobian block but not the
reporting block; reporting directional-Jacobian evidence remains a
promotion-receipt item rather than motivating a new runtime seam.
