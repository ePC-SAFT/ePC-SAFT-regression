# Pure-Saturation Regression

The package fits segment count, segment diameter, and dispersion energy for
one admitted pure component per call. It reproduces the accepted methane
workflow and adds a local ethane candidate. The installed `epcsaft` provider
owns Helmholtz energy and nonlinear derivatives.

## Source records

The methane table contains NIST Chemistry WebBook SRD 69 saturation properties
from 100 through 180 K in 10 K increments. The fit uses 110, 130, 150, and
170 K. The other five rows remain descriptive held-out evidence. The retained
lab file has SHA-256
`a5e16df3bf8ec78483fc340782cddc89ab8b284a9f6dfaecd6cda3ffde579227`.
The package normalizes CRLF to LF and has SHA-256
`dec64d5a6cac414a4a92393a0d728fa27c02135c6a159d0d1881d7b6dde6d26c`.

The ethane table contains NIST SRD 69 saturation properties from 100 through
280 K in 20 K increments. The fit uses 140, 180, 220, and 260 K. The held-out
set uses 120, 160, 200, and 240 K. The package reports 100 and 280 K as domain
stress. The retained CRLF file has SHA-256
`ed09b8781acfb7025ca505878b884f6353ddd9f3f4bd7aae2e6df88bbe847a67`.
The LF-normalized package file has SHA-256
`b01333e827933c0a7148672c8ae3eef78393320c0d18f2c4d5a0fc40d9bef6b2`.

Both tables keep temperature in K, saturation pressure in Pa, and saturated
liquid mass density in kg/m3. The records retain the exact NIST query, decimal
strings, retrieval date, and use basis.

## Provider contract

Each phase evaluation fixes `n = 1 mol` and sends

```text
(n_mol, V_m3, m, sigma_angstrom, epsilon_over_k_kelvin)
```

to `evaluate_pure_phase_parameters`. The provider returns
`Phi = A/(RT n_ref)`, its five-coordinate gradient and Hessian, pressure,
chemical potential divided by `RT`, and its immutable parameter fingerprint.
Regression evaluates

```text
dP/dq_j       = -R T Phi_(V,q_j)
d(mu/RT)/dq_j = Phi_(n,q_j)
```

and owns no EOS equation.

## Component specifications

Methane keeps the accepted start and scale:

```text
start = (1.08, 3.555744 angstrom, 157.5315 K)
scale = (0.1, 0.1 angstrom, 10 K)
```

Ethane starts from the Gross 2001 provider parameters:

```text
start = (1.6069, 3.5206 angstrom, 191.42 K)
scale = (0.1, 0.1 angstrom, 10 K)
```

Both specifications use bounds `[0.5, 3.5]`, `[2, 5] angstrom`, and
`[50, 400] K`. Methane fixes molar mass at `0.016043 kg/mol`; ethane fixes it
at `0.030070 kg/mol`. Each specification also fixes its source, dataset,
training partition, provider fingerprint, volume bounds, pressure bounds,
solver controls, and confirmation gate.

The ethane vapor-volume upper bound is `100 m3` and its reporting-pressure
lower bound is `1 Pa`. Those values contain the audited low-temperature rows.
They do not relax residual closure or topology thresholds.

## Fit formulation

The physical parameters use affine dimensionless coordinates:

```text
p = start + scale*z
```

Each training row carries liquid-like and vapor-like volumes. Log transforms
make both positive. Disjoint volume bounds enforce liquid volume below vapor
volume. Observed liquid density and ideal-gas volume provide starts.

Each row contributes liquid pressure, vapor pressure, liquid-vapor chemical
potential equality, and liquid mass-density residuals. Pressure uses observed
pressure as scale. Chemical potential uses scale 1, and density uses observed
density. Each residual receives weight 0.25 for numerical row normalization.
The weights make no uncertainty claim.

The native cost assembles the `16 x 11` Jacobian from the provider Hessian and
the affine, logarithmic, density, and scaling chain rules. A parameterized
test checks methane and ethane directional products against centered residual
differences. Production exposes no numerical derivative option.

## Acceptance layers

Solver convergence requires Ceres `CONVERGENCE`, a usable finite solution,
complete Jacobian columns, full parameter-column rank, parameter movement
inside the declared bounds, and no provider error. A confirmation solve
perturbs the volume starts and requires scaled parameter agreement below
`1e-5` and relative cost agreement below `1e-8`.

Physical acceptance requires positive stability slopes, ordered distinct
phases, relative phase-volume separation above `1e-3`, and reporting closure
on training and held-out rows. Pressure closure uses `1e-8` after scaling by
observed pressure. Chemical-potential closure uses `1e-8`.

Stress rows retain their termination, closure, topology, predictions, and
failure reasons. They do not enter solver or physical acceptance and do not
receive aggregate errors unless the reporting state passes the physical
gates. In the local ethane candidate, 280 K passes. The 100 K solve terminates
at the `1 Pa` lower pressure bound and fails pressure and chemical-potential
closure, so the receipt records no valid 100 K prediction.

The result reports full and parameter-column singular values, ranks, and
condition numbers. Those local diagnostics provide no global identifiability
or uncertainty claim. No source-row error has an approved pass threshold.

## Failure boundaries

The package rejects source, component, unit, row-order, partition,
specification, and provider-fingerprint mismatches before Ceres starts. Native
callbacks retain compact provider, topology, derivative, solve, and reporting
failures.

The compiled problem validates the full source and specification identity and
returns it to Python. Runtime code claims no provider commit or wheel digest.
The candidate runner binds and byte-checks the installed wheels before import.

The private native test surface exposes the training residual and Jacobian
block. It does not expose the reporting Jacobian, and this slice adds no test-
only runtime seam.

The package owns no binary interaction, association, electrolyte, reactive,
generic-family, persistence, or provider-catalog capability.
