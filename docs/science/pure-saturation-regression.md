# Pure-Saturation Regression

The package fits segment count, segment diameter, and dispersion energy for
one closed pure-component record per call. Methane and ethane are accepted
reproducible workflows. Propane Checkpoint A is an authority-neutral local
candidate. The installed `epcsaft` provider owns Helmholtz energy and nonlinear
derivatives.

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

The propane table is the exact 24-row target from the approved Glos,
Kleinrahm, and Wagner 2004 direct-experimental packet at Validation commit
`7e51590757f1cb85f51df98e9fe1f88cd4255a88`, tree
`05af9e948c786ddfcf43dba701970f1cbb6435a2`. The packaged target SHA-256 is
`ccd1cfa15ec44432b06cbf22316d168c61b282631c9b1e1591e497b8d48b5676`;
the governing packet YAML is
`ba31448989f565d05d63908076e836977780aa87199f208310e9b80b03f64697`;
and the 63-row source receipt is
`ed5eb703ccd3e6bb4c4cfa82ecd58c58f9da0c93ab07a204dee94d8b0ae8d081`.
Training uses 150, 210, 270, and 330 K. The 18 remaining interior rows are
held out, while 110 and 340 K are stress rows. Pointwise pressure, liquid-
density, and available vapor-density uncertainties are retained on the source
observations. They do not weight residuals or define model-acceptance cutoffs.

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

Propane starts from the provider's source-backed Gross 2001 propane bundle:

```text
start = (2.0020, 3.6184 angstrom, 208.11 K)
scale = (0.1, 0.1 angstrom, 10 K)
provider fingerprint = sha256:9bfbc8d7789e51609945e61dbdf7a020decc8f9e31b408b0977724c7cb3e1551
```

Both specifications use bounds `[0.5, 3.5]`, `[2, 5] angstrom`, and
`[50, 400] K`. Methane fixes molar mass at `0.016043 kg/mol`; ethane fixes it
at `0.030070 kg/mol`; propane fixes it at `0.044096 kg/mol`. Each
specification also fixes its source, dataset, training partition, provider
fingerprint, volume bounds, pressure bounds, solver controls, and confirmation
gate.

The ethane vapor-volume upper bound is `100 m3` and its reporting-pressure
lower bound is `1 Pa`. Those values contain the audited low-temperature rows.
They do not relax residual closure or topology thresholds.

Propane uses liquid-volume bounds `[2e-5, 1.2e-4] m3`, vapor-volume bounds
`[1.5e-4, 2000] m3`, and reporting-pressure bounds `[0.1, 1e7] Pa`. These
domains contain every source-derived liquid-volume and ideal-gas vapor-volume
start, including the 0.6 Pa stress row. They are solver domains, not accuracy
tolerances. The propane iteration budget is 5000 because the unchanged Ceres
tolerances converge at iterations 1090 and 1165 for the primary and perturbed
starts; the accepted methane and ethane budgets remain 500 and their numerical
results are unchanged.

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
test checks methane, ethane, and propane directional products against centered
residual differences. Production exposes no numerical derivative option.

## Acceptance layers

Solver convergence requires Ceres `CONVERGENCE`, a usable finite solution,
complete Jacobian columns, full parameter-column rank, parameter movement
inside the declared bounds, and no provider error. A confirmation solve
perturbs the volume starts and requires scaled parameter agreement below
`1e-5` and relative cost agreement below `1e-8`. Cost agreement is the
symmetric relative difference `abs(primary - confirmation) /
max(abs(primary), abs(confirmation), numeric_limits<double>::min())`, so the
gate remains relative when both fitted costs are below one.

Physical acceptance requires positive stability slopes, ordered distinct
phases, relative phase-volume separation above `1e-3`, and reporting closure
on training and held-out rows. Pressure closure uses `1e-8` after scaling by
observed pressure. Chemical-potential closure uses `1e-8`.

Stress rows retain their termination, closure, topology, predictions, and
failure reasons. They do not enter solver or physical acceptance and do not
receive aggregate errors unless the reporting state passes the physical
gates. In the accepted ethane workflow, 280 K passes. The 100 K solve terminates
at the `1 Pa` lower pressure bound and fails pressure and chemical-potential
closure, so the receipt records no valid 100 K prediction.

For propane, the primary and perturbed training solves converge with full rank
and agree inside the unchanged confirmation gates. The 110 K stress solve
reaches the vapor-volume ceiling and has no acceptance effect. The 120 K
held-out reporting solve converges and is usable, but its liquid-pressure
closure is approximately `1.0540e-7 Pa`, or `3.29e-8` after division by the
observed `3.2 Pa`; this exceeds the frozen `1e-8` numerical closure gate.
Tightening a diagnostic-only reporting step tolerance did not improve the
provider-returned pressure plateau and was not retained. The propane result is
therefore solver-converged and numerically confirmed but not physically valid.
Predictive status remains `NOT_ADJUDICATED_NO_APPROVED_HELD_OUT_CUTOFF`; neither
the source uncertainty nor this numerical failure is converted into a model-
accuracy decision.

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
