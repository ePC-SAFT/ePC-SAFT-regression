# Pure-Methane Saturation Regression

This package owns one candidate regression formulation. It fits methane
segment count, segment diameter, and dispersion energy from four NIST
saturation rows. The installed `epcsaft` provider owns Helmholtz energy and its
nonlinear derivatives.

## Source records

The retained table contains NIST Chemistry WebBook SRD 69 methane saturation
properties at 100–180 K in 10 K increments. The package keeps temperature in
K, saturation pressure in Pa, and saturated-liquid mass density in kg/m3. It
uses 110, 130, 150, and 170 K for fitting and reports all nine rows.

The source query identifies methane by CAS 74-82-8 and requests eight-digit
output. The retained lab artifact has SHA-256
`a5e16df3bf8ec78483fc340782cddc89ab8b284a9f6dfaecd6cda3ffde579227`.
The package changes CRLF line endings to LF and records the packaged hash in
`epcsaft_regression.records`. Decimal strings and fields do not change. NIST
SRD redistribution and use terms apply.

## Provider contract

Each phase evaluation fixes `n = 1 mol` and sends

```text
(n_mol, V_m3, m, sigma_angstrom, epsilon_over_k_kelvin)
```

to `evaluate_pure_phase_parameters`. The provider returns
`Phi = A/(RT n_ref)`, its five-coordinate gradient and 5×5 Hessian, pressure,
chemical potential divided by `RT`, and its immutable source-parameter
fingerprint. Regression uses

```text
dP/dq_j       = -R T Phi_(V,q_j)
d(mu/RT)/dq_j = Phi_(n,q_j)
```

and does not evaluate, copy, or differentiate any EOS equation.

## Fit formulation

The physical parameters use affine dimensionless coordinates:

```text
p = start + scale*z
start = (1.08, 3.555744 angstrom, 157.5315 K)
scale = (0.1, 0.1 angstrom, 10 K)
bounds = ([0.5, 3.5], [2, 5] angstrom, [50, 400] K)
```

The start is retained from the M5 plan as the deliberate displacement
`(1.0*1.08, 3.7039*0.96, 150.03*1.05)`. The fixed methane molar mass
`0.016043 kg/mol` is also an admitted formulation constant from that plan.

Each training row carries liquid-like and vapor-like volumes. Log transforms
make both positive. Disjoint volume bounds enforce liquid volume below vapor
volume. The observed liquid density and ideal-gas volume provide starts; they
do not remove either volume from the problem.

The four raw residuals per row are liquid pressure minus observed pressure,
vapor pressure minus observed pressure, liquid-vapor chemical-potential
difference, and calculated minus observed liquid mass density. Pressure uses
the observed pressure as scale. Chemical potential uses scale 1, and density
uses the observed density. Each residual receives weight 0.25 for equal row
normalization. The weights make no uncertainty claim.

The native cost assembles the 16×11 Jacobian from the provider Hessian and the
exact affine, logarithmic, and density chain rules. Test code checks a fixed
directional product with centered residual differences. Production exposes no
numerical derivative option.

## Acceptance layers

Solver convergence requires Ceres `CONVERGENCE`, a usable finite solution,
complete Jacobian columns, and respected bounds. A second solve perturbs all
row-volume starts and checks scaled parameter and cost agreement. Physical
validity requires positive `dP/d(rho_m)`, ordered distinct volumes, finite
diagnostics, relative phase-volume separation above `1e-3`, and successful
saturation reporting solves at all nine source temperatures. The reporting
result reapplies these gates after Ceres terminates; a rejected final row keeps
its row-specific failure reason.

The confirmation volume-start multipliers are `1.01` for liquid and `0.98` for
vapor; its scaled-parameter and relative-cost limits are `1e-5` and `1e-8`.
Ceres uses `DENSE_QR`, one thread, silent logging, 500 iterations maximum, and
three `1e-10` termination tolerances.

Each reporting solve uses `P=P_observed*exp(u_P)` bounded to `[1e3,1e7] Pa`.
Pressure closures are scaled by observed pressure and must be at most `1e-8`;
the dimensionless chemical-potential closure must also be at most `1e-8`.
These are numerical equation-closure gates, not measurement uncertainties.

The result reports full and parameter-column singular values, ranks and
condition numbers. These local diagnostics do not prove global
identifiability. Reporting errors have no pass threshold in this slice;
installed-artifact validation owns predictive admission.

## Failure boundaries

The package rejects changed source identity or units, duplicate or reordered
rows, nonfinite values, unsupported species, capsule ABI mismatch, missing
parameterized transport, provider evaluation failure, unstable phases,
coalescence, topology loss, incomplete derivatives, nonconverged solves, and
seed-dependent results. Provider, topology, and post-solve failures for valid
inputs return compact native evidence in `MethaneFitResult.failure_reasons`;
contract parsing errors may still raise before a solve outcome exists.

The compiled native problem independently validates and retains all exact row
and source identities, units, transforms, bounds, solver controls, provider
SDK ABI/table contract, and provider source fingerprint. Native diagnostics
return their own row IDs, and the full compiled identity must round-trip
unchanged. The runtime does not claim a provider commit or wheel digest. The
isolated candidate runner hashes and byte-binds those exact artifacts before
import, and the canonical receipt retains that artifact evidence.

Training exact-Jacobian directional coverage uses the existing private native
test surface. That surface does not expose the reporting block Jacobian; its
directional check remains a promotion-receipt evidence item, avoiding a new
production seam solely for testing.

The package does not own binary interactions, association, electrolytes,
reactions, generic targets, persistence, or provider parameter catalogs.
