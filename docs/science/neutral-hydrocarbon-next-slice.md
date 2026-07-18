# Neutral-Hydrocarbon Next-Slice Contract

Status: binary preflight blocked by frozen pressure closure; authority effect
none.

This record freezes the smallest regression-owned contract for two ordered
checkpoints. It adds no runtime fit, public export, provider equation,
provider-catalog parameter, promotion, or predictive claim.

## Checkpoint order and readiness

Checkpoint A extends the existing pure-saturation engine to propane. Runtime
work remains at
`READY_WAITING_PROPANE_PRIMARY_SOURCE_PACKET` until validation commits its new
direct-experimental primary-source packet. Validation commit `267f853` retains
NIST WebBook reference-EOS calculations, so it does not meet this gate. No lab
file can substitute for experimental source evidence.

Checkpoint B fits one constant methane/ethane `kij`. The manager-retained
provider subject is commit
`45d5764f61729d387100348a8ff91792f6e0a395`, tree
`271d4848faf73afd4cf0683efe5c855053df7d01`, wheel SHA-256
`95d2292b052ab74657931f2dec97c3ea4160d9b17812956515c7195a853e6c5b`,
and installed public-header SHA-256
`6f3a186bf5359f32449a31c544e8b8525be6c594804151d3e74a74e411ded8f4`.
The exact installed artifact passes the derivative and numerical feasibility
checks below, but runtime work is `BLOCKED_FROZEN_PRESSURE_CLOSURE`: the
approved fixed-composition formulation misses its `1e-8` physical closure
gate. The preflight is a non-installed evidence target, not a package API or
runtime implementation.

## Frozen source and model identity

The binary target uses all 17 rows in validation's audited May 2015 packet:

- validation source commit:
  `73a37f5935e919a34d1e4fa3af285951d6fac8e7`;
- source CSV SHA-256:
  `5cd1e74925a3c6504f5106dcf911f2cae2d6e99a5133fccc20454d8991bdbc7f`;
- source packet SHA-256:
  `d43433e93b354e01f96d330c760818a24b775026461ce795e45774cfb11ac94e`;
- component order: methane, ethane; and
- expected base-provider fingerprint:
  `sha256:307fcb28d535b94782f3e3caf4012c0c8c0dc87ee4239d6c316de56553543286`.

The clean Gross--Sadowski pair record fixes the source value `kij = 0` and
identifies the mixing rule used by the provider as

```text
epsilon_ij = sqrt(epsilon_i * epsilon_j) * (1 - kij).
```

The pair-record file SHA-256 is
`747e8281c7a1e4240ee4badbc0bedd047521fb303726699b38c99fccf7f74c2a`.
The regression start is therefore exactly zero. The approved engineering
search interval is `[-0.15, 0.10]`; it is a bounded neighborhood of the source
value and keeps the cross-dispersion energy positive. It is not attributed to
the paper. The transformed parameter is

```text
kij = 0 + 0.01 * z_kij.
```

The fit is rejected if `kij` is nonfinite or within `1e-8` of either bound.
Bounds cannot be widened after observing the fit.

## Variables and scaling

For row `r`, measured `T_r`, `P_r`, liquid methane fraction `x_r`, and vapor
methane fraction `y_r` remain fixed. Each phase contains one mole:

```text
n_L = (x_r, 1 - x_r) mol
n_V = (y_r, 1 - y_r) mol
```

Each row lifts two positive volumes:

```text
V_L = 6.5e-5 m3 * exp(u_L)
V_V = (R * T_r / P_r) * exp(u_V)
```

The liquid scale lies inside the `5.54e-5` to `8.56e-5 m3` range retained by
the current clean, authority-neutral May campaign for its 16 evaluated local
states. It is initialization evidence only. The vapor scale is the
dimensionally exact ideal-gas one-mole volume. Frozen bounds are
`[2e-5, 1e-4] m3` for liquid and `[1e-4, 1e-2] m3` for vapor. They contain the
current source-state initialization evidence and enforce phase ordering
without density roots.

The complete variable order is one `z_kij`, followed by `(u_L, u_V)` in the
17 immutable source-row order. The problem therefore has 35 dimensionless
variables.

## Residuals

The provider owns

```text
Phi = A / (R T n_ref),  n_ref = 1 mol
P = -R T n_ref * Phi_V
mu_i / (R T) = Phi_ni
```

For each row, regression assembles four residuals with common numerical weight
`w = 0.25`:

```text
r_LP  = w * (P_L - P_r) / P_r
r_VP  = w * (P_V - P_r) / P_r
r_mu1 = w * (Phi_n1_L - Phi_n1_V)
r_mu2 = w * (Phi_n2_L - Phi_n2_V)
```

All 68 residuals are dimensionless. The two pressure residuals use the one
observed coexistence pressure to define each lifted phase state; they are not
treated as two independent measurements. The weights make no uncertainty
claim. Measured compositions are inputs, not lifted data-reconciliation
variables.

## Exact Jacobian contract

The provider callback must return `Phi`, its gradient, and its row-major
symmetric Hessian in exact coordinate order

```text
(n_methane_mol, n_ethane_mol, V_m3, kij).
```

For phase `a` and component `i`, the consumed derivatives are

```text
d(r_aP)/d(u_a)   = w/P_r * (-R*T*n_ref*Phi_VV_a) * V_a
d(r_aP)/d(z_kij) = w/P_r * (-R*T*n_ref*Phi_Vk_a) * 0.01

d(r_mui)/d(u_L)   =  w * Phi_niV_L * V_L
d(r_mui)/d(u_V)   = -w * Phi_niV_V * V_V
d(r_mui)/d(z_kij) =  w * (Phi_nik_L - Phi_nik_V) * 0.01
```

All other row-to-volume entries are zero. This produces the exact `68 x 35`
Ceres residual Jacobian from a second-order CppAD provider tensor. No third
derivative, density root, implicit branch derivative, equilibrium dependency,
or numerical production derivative is consumed.

## Future result and acceptance diagnostics

These are the documented requirements for the later binary implementation;
no forward result or diagnostic types are packaged by the current preparation.

The private immutable result contract retains exact source/provider hashes;
separate solver, numerical, and physical statuses; termination and usability;
initial/final costs; iteration count; start, final, movement, bounds, and
active-bound state for `kij`; every source row and its raw/scaled residuals;
phase volumes, pressures, chemical potentials, stability slopes, and rejection
reasons; full and projected Jacobian diagnostics; confirmation outcomes; and
the predictive status
`NOT_ADJUDICATED_NO_APPROVED_HELD_OUT_CUTOFF`.

Solver acceptance requires finite complete Jacobian columns, Ceres
`CONVERGENCE`, a usable finite solution, cost reduction from the source start,
full scaled-Jacobian rank 35, projected-parameter rank 1, finite nonzero
parameter movement, and a finite non-bound fitted `kij`. Rank uses

```text
sigma_max * max(68, 35) * epsilon_double.
```

The projected parameter column is

```text
j_k_perp = (I - J_V * pinv(J_V)) * j_k.
```

Its one nonzero singular value and rank are retained; no covariance,
uncertainty, or global-identifiability interpretation is allowed.

Numerical confirmation perturbs all liquid starts by `1.01`, all vapor starts
by `0.98`, and repeats from `kij = -0.05` and `+0.05`. It requires maximum
scaled-`kij` agreement `<= 1e-5` and symmetric relative-cost agreement
`<= 1e-8`.

Physical status requires exact source and provider identity, positive amounts
and volumes, provider-domain success, ordered separated phases, stable phase
slopes, and scaled pressure closure `<= 1e-8`. The fitted chemical-potential
residuals remain reported calibration errors because no accuracy cutoff was
approved.

## Installed-artifact preflight result

The non-installed preflight consumed the exact retained wheel and installed
header, the audited CSV above, and no sibling provider source. Its deterministic
receipt is `evidence/binary-kij-preflight.json` (subject SHA-256
`d51c9f0713b6a7355be719b6843e4459f41d46d16d973668d694715d36b63676`, file
SHA-256
`27c7eab578290a791f47a4aabd16d4a4b95949d13a6b6b022d482c5130309f28`).

Provider value/gradient/Hessian directional checks pass, including the active
`kij` coordinate. The exact scaled residual Jacobian directional check has
maximum absolute error `3.2377076436551633e-9` and scaled error
`0.012687348470964473`. The primary solve converges in seven iterations to
`kij = -0.008430339623289762`, away from both declared bounds. Its complete
scaled Jacobian has rank 35, condition number `1345.943122509698`, and rank
threshold `3.814411826697633e-13`; the projected parameter singular value is
`0.01886370066465226` and its rank is 1.

The starts at `-0.05` and `+0.05` converge to `-0.008430317065561237` and
`-0.00843032253373879`. Their scaled parameter deltas from the primary result
are `2.255772852538529e-6` and `1.7089550972265721e-6`; symmetric relative
cost deltas are `6.045202596244109e-11` and `7.31892665322379e-11`.

Solver convergence and numerical confirmation therefore pass. Physical
acceptance does not: the maximum relative pressure closure is
`0.022691483584976503`, not `<= 1e-8`. The worst row is
`may2015-ch4-c2h6-002`, with liquid and vapor relative pressure closures
`0.003910995749513778` and `-0.022691483584976503`, methane and ethane
`mu/RT` equality residuals `0.014205213221523572` and
`-0.06950599355260767`. No tolerance, residual weight, variable, or target was
changed after observing this result.

## Explicit exclusions

This preparation does not admit predicted `x/y`, held-out binary evidence,
uncertainty, covariance, temperature-dependent interactions, other binary
pairs, association, electrolytes, equilibrium, density closure, global
identifiability, parameter persistence, a generic target registry, promotion,
publication, or authority transfer.
