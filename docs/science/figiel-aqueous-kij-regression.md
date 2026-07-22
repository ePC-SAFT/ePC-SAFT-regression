# Figiel Aqueous Current-Catalog Interaction Regression

Status: installed derivative/rank/multistart preflight complete; runtime is
blocked by two frozen gates. Exact callback columns 1, 3, 4, and 5 fail the
predeclared directional check. The rank-11 bounded solve converges consistently
from all three starts, but its largest published-parameter difference is `1.8`
and five parameters are on bounds. No production Ceres fit, Regression wheel,
Validation campaign, parameter admission, prediction, or authority transfer
exists for this slice.

This document is the sole Regression design owner for the next Figiel parameter
family after the five-ion Born-diameter tracer. It implements the staged intent
of Migration D-023 without creating a generic parameter registry or combining
later organic-solvent and expanded-ion work.

## Claim and scope

The bounded claim is that one source-bound Ceres fit can recover the 11
current-catalog aqueous interaction parameters printed as fitted in Figiel,
Yu, and Held (2025) Tables 4 and 5 while reproducing the same 164 aqueous
mean-ionic-activity-coefficient observations used as training evidence.

The active vector, in immutable order, is

| Column | Interaction | Published `k_ij` | Source cell |
| ---: | --- | ---: | --- |
| 0 | water--Li+ | -0.4 | Table 5, Li+ row, water column |
| 1 | water--Na+ | -0.3 | Table 5, Na+ row, water column |
| 2 | water--K+ | -0.1 | Table 5, K+ row, water column |
| 3 | water--Cl- | -0.3 | Table 5, Cl- row, water column |
| 4 | water--Br- | -0.3 | Table 5, Br- row, water column |
| 5 | Li+--Cl- | 0.8 | Table 4, Cl- row, Li+ column |
| 6 | Na+--Cl- | 0.8 | Table 4, Cl- row, Na+ column |
| 7 | K+--Cl- | 0 | Table 4, Cl- row, K+ column |
| 8 | Li+--Br- | 0.5 | Table 4, Br- row, Li+ column |
| 9 | Na+--Br- | 0.65 | Table 4, Br- row, Na+ column |
| 10 | K+--Br- | -0.35 | Table 4, Br- row, K+ column |

The explicit fitted zero is a parameter target. Blank Table 4 or 5 cells are
not zeros and are excluded. Water--methanol and water--ethanol interactions are
inherited rather than fitted in Figiel and are excluded.

The slice does not include H+, I-, sulfate, vanadium ions, methanol, ethanol,
solvent solvation factors, dielectric suppression, Born diameters, neutral-
solvent parameters, density or osmotic targets, uncertainty, prediction,
catalog persistence, or a simultaneous Tables 2--5 solve.

This is a conditional `k_ij` family recovery. Every non-`k_ij` model input is
fixed to the immutable published Figiel Provider catalog: Table 2 water
parameters including `f_water = 1.5`, Table 3 ion parameters and Born
diameters, and the equation-11 dielectric coefficient `7.01`. The locally
recovered Born-diameter candidate is not a Provider-catalog artifact and is not
substituted through an overlay. This checkpoint therefore proves the active
aqueous interaction family only; later stage-C recovery and the final joint
tuple replay must test cross-family coupling separately.

## Sources and existing forward evidence

The primary paper is Figiel, Yu, and Held, *Industrial & Engineering Chemistry
Research* 64 (2025) 9406--9418, DOI `10.1021/acs.iecr.5c00475`. The retained
read-only Markdown is 71,826 bytes with SHA-256
`ce80533925a91bc59d8d0d8056113c40611ca26c2edf04aced76986d50bd4bae`.
Its Model Parameters section states that aqueous ion--water and ion--ion
interactions were adjusted to experimental MIAC data at 298.15 K; Tables 4 and
5 provide the values above.

The approved Validation packet is commit
`8944d34f7002cda1bb8760e606cc1f11696f58cd`, tree
`6c8fd350dcd6bfdd7be1918f73fd33a23e2070dd`. Its relevant identities are:

- target ledger SHA-256
  `f405a3e48d21cd979a8dd480d5f8cb3be40754f5d6babf368b505b5f305607f0`;
- parameter-provenance SHA-256
  `932e8baa90fcefbaa8c3a8730cdeadd83a4c01f0a3b109f4e4cd0319aee9312b`;
- metadata SHA-256
  `8ea06c6ca5452d01448a03f9a76cf7d0c35bb99c9abe23ccb1729d56c71d468f`;
  and
- Hamer--Wu 1972 aqueous alkali-halide CSV SHA-256
  `2f63e13f06a5b0f4e8bca2980b6a8d9d7fb0f839153c43e3a71952daf9796595`.

The Hamer--Wu packet contains 164 positive molality-scale MIAC rows for LiCl,
NaCl, KCl, LiBr, NaBr, and KBr at 298.15 K and 1 bar. All rows are training
observations. The paper reports no rowwise uncertainty and does not disclose a
reproducible objective, weights, parameter bounds, or starts; this design makes
those Regression choices explicit and does not attribute them to the paper.

Validation's existing published-tuple forward campaign is useful baseline
evidence, not parameter-recovery evidence. Its JSON SHA-256 is
`411d404a4f4f7cec28be065af69d7e62d2e9014fc4d3ab545519f1bb63d1aa44`.
With the published tuple it reports pooled MIAC RMSE `0.13407022131157917` and
passes its bounded physical-reproduction campaign. It contains no active
parameter derivatives, Ceres solve, rank result, multistart result, or recovered
parameter vector.

## Residual and variables

For observation `q`, with positive observed molality-scale mean ionic activity
coefficient `gamma_q_obs`, Regression owns the dimensionless residual

```text
r_q(k) = ln(gamma_q_model(k) / gamma_q_obs).
```

All 164 residuals have equal weight because the source packet has no pointwise
uncertainty. The logarithmic ratio treats equal multiplicative deviations
symmetrically and does not require a dimensional scale.

The Ceres problem has 11 variables and 164 residuals. Each parameter is
dimensionless with transform `k_j = z_j`, scale 1, and bounds `[-1, 1]`. The
bound is a source-enclosing engineering interval for all Figiel fitted
interactions and keeps the dispersion combining factor `1-k_ij` nonnegative.
The current 11 published targets occupy `[-0.4, 0.8]` and are interior.

The deterministic starts are:

1. all `0` (primary);
2. all `-0.5`; and
3. all `+0.5`.

The published tuple is comparison evidence only and is not a solver seed.
Ceres uses `DENSE_QR`, one thread, silent logging, at most 500 iterations, and
function, gradient, and parameter tolerances of `1e-10`.

## Exact Provider contract

Provider implementation `8ae37dbc4dc61a4ee1109bb1cf8e26470e457975`,
merged as `37bcece2b608421434df044be4c5e2b0c67b946d`, supplies the required
installed callback. The retained wheel SHA-256 is
`d984d5b68d97d5a40b8d4d729d8ceff5aa0e274222ca548d6fdf82162869537f` and
the installed public-header SHA-256 is
`01568808f48c8cf0cd5fd0eb0b3d038349a319251eb8a9f9053b028ed35e5a36`.

The minimum new Provider contract is one model-bound aqueous-alkali-halide
evaluation. For the selected ordered model `(water, cation, anion)`, fixed
`T`, fixed `P`, and one formula-unit molality, it accepts exactly three finite
trial interactions in this order:

```text
(k_water_cation, k_water_anion, k_cation_anion).
```

It returns:

- `ln(gamma_pm^m)`;
- its three exact total fixed-`T,P` derivatives with respect to those active
  interactions;
- the terminal reference molality and convergence diagnostic;
- the parameter fingerprint and component order; and
- structured unsupported-model, domain, numerical, and callback status.

Provider retains the reference sequence, density closure, EOS, combining rule,
and CppAD ownership. The total first parameter derivative requires mixed
state/parameter second derivatives of the Provider Helmholtz owner and the
first implicit density sensitivity. It requires no third derivatives. Ceres
consumes only the returned first total derivatives. Regression must not copy
the reference sequence, solve density independently, modify parameter records,
or use production finite differences.

For a row belonging to salt `(c,a)`, Regression maps the three returned
derivatives into the global water--cation, water--anion, and cation--anion
columns. Every other entry in that row is exactly zero. The assembled Jacobian
therefore has shape `164 x 11` with exactly three structurally active columns
per row.

## Preflight and identifiability

Before runtime implementation, an installed-artifact preflight evaluates the
complete residual and Jacobian at one nonpublished interior trial vector. It
checks every global column against centered callback-value differences at
steps `1e-4` and `5e-5`. For each column the infinity-norm difference must meet

```text
abs(J_exact - J_h/2)
  <= max(1e-8, 20*abs(J_h - J_h/2), 2e-8*abs(J_exact)).
```

Finite differences are preflight evidence only. They are not a production
backend.

SVD uses rank threshold

```text
s_max * max(164, 11) * epsilon_binary64 * 100.
```

The preflight must show 11 complete nonzero columns and rank 11. It retains all
singular values, the condition number, normalized sensitivity-column `J^T J`
correlations, and the least-sensitive right singular vector. These are local
sensitivity diagnostics, not a covariance or uncertainty estimate. If rank is
below 11, this exact joint fit stops. Regression does not add regularization,
priors, parameter tying, or fixed published values to manufacture rank.

The installed-artifact preflight uses the nonpublished all-`0.2` trial. It has
11 complete nonzero columns, rank 11 at threshold
`3.396827102277414e-11`, condition number `198577.925594654`, terminal
reference molality `1e-12`, maximum reference convergence diagnostic
`1.0628418010583118e-5`, and maximum derivative convergence diagnostic
`8.281553220967908e-11`. Columns 0, 2, and 6--10 pass the two-step directional
check; columns 1, 3, 4, and 5 fail it. This is an exact-derivative evidence
failure, even though every column is finite, nonzero, and the global matrix is
full rank.

The bounded nonlinear preflight then exercises the complete frozen objective
without creating package runtime. All-zero, all-`-0.5`, and all-`+0.5` starts
converge to the same rank-11 vector within
`1.0952261320085199e-10`. The primary cost is
`0.25101017330977715`, below the published tuple's
`0.4200114112124652`, but the fitted vector differs from the printed tuple by
as much as `1.8`; columns 3, 5, 6, 7, and 8 are on the lower bound. It therefore
fails both the `0.05` recovery gate and the non-bound gate. This is stronger
than the earlier published-point gradient observation, but it does not prove
that every possible source-backed objective would miss the table.

Figiel does not disclose the objective, weights, bounds, starts, or fitting
sequence used to obtain the printed values. Regression must not tune those
choices post hoc, seed from the printed tuple, or add a prior merely to
reproduce the answer. Exact compact evidence is retained in
`evidence/figiel-aqueous-kij-published-tuple-preflight.json`, SHA-256
`12cb7205e988316f9c61560bc82e012b48d681fa9107e6da538b142fe695078c`.

## Result and acceptance semantics

One closed result retains the source and artifact identities, all 11 fitted
parameters and published differences, 164 row diagnostics, three starts,
objective values, gradients, rank/conditioning/sensitivity diagnostics,
active-bound status, Provider diagnostics, and ordered failure reasons.

The statuses remain separate:

1. `solver_converged`: primary Ceres termination, usability, finiteness, cost
   reduction, callback, and bound conditions pass.
2. `numerically_converged`: all three starts satisfy the solver gates, rank is
   11, and their fitted vectors agree within `1e-5` in infinity norm.
3. `physically_valid`: every fixed-pressure reference and target state is
   Provider-usable, every returned `ln(gamma_pm^m)` is finite, and its implied
   activity coefficient is positive.
4. `workflow_valid`: physical validity plus source hashes, exact 164-row
   membership, units, salt and component order, state, fixed non-`k_ij`
   catalog basis, Provider artifact, fingerprints, and reference diagnostics
   match the frozen contract.
5. `scientifically_valid`: workflow validity plus the parameter-proximity and
   observable-reproduction gates below pass.
6. `predictive_status` is
   `NOT_ADJUDICATED_NO_APPROVED_HELD_OUT_CUTOFF` because every row is training
   evidence.

The user-approved parameter-recovery gate is

```text
max_j abs(k_j_fit - k_j_published) <= 0.05.
```

For `u_ij = sqrt(u_i*u_j)*(1-k_ij)`, this limits the recovered change to five
percent of the geometric-mean cross-dispersion energy scale. It is an
engineering recovery criterion, not source uncertainty or a claim that the
printed parameters are unique.

The fitted log-residual cost must not exceed the cost of the published tuple,
allowing only `C_fit <= C_published*(1 + 1e-10) + 1e-14` for binary64
comparison. The installed-artifact Validation replay must also retain the
existing forward-reproduction checks:

- pooled MIAC RMSE no greater than `0.17`;
- each salt RMSE no greater than `0.35`;
- each salt maximum absolute MIAC error no greater than `1.25`; and
- the first predicted MIAC for each salt below `0.98`.

These are in-sample physical-reproduction checks, not experimental uncertainty
or predictive cutoffs. A converged fit that fails rank, multistart, parameter
proximity, or observable reproduction remains scientifically invalid.

## Ownership and execution order

The eventual implementation reuses the existing Python package, immutable
record owner, workflow/result owner, native module, Ceres dependency, and one
CMake target. It may add one family-specific native cost owner; it must not add
a second module or target, generic target registry, parameter overlay, mutable
catalog, compatibility layer, copied Provider equation/reference logic, or
Equilibrium dependency.

The executable order is:

1. Preserve the exact Provider wheel/header and complete blocked preflight.
2. Provider diagnoses exact derivative columns 1, 3, 4, and 5 at the retained
   all-`0.2` trial without changing Regression's evidence criterion.
3. After derivative correction or confirmation, resolve the scientific claim:
   keep the equal-weight objective and report failed table recovery, or supply
   a source-backed alternative. Do not tune weights, bounds, starts, or
   tolerances to the printed answer.
4. Use bounded independent subagent review on the exact revised contract.
5. Only if every frozen gate passes, implement the closed Ceres fit and retain one commit-bound wheel
   plus candidate evidence.
6. Regression then authors the bounded public-installed-artifact campaign in
   Validation under Migration's serialized writer protocol.
7. A distinct independent review decides admission. No package-authored result
   self-promotes or writes fitted values into the Provider catalog.

The present state is
`BLOCKED_FROZEN_DERIVATIVE_AND_PARAMETER_RECOVERY_GATES`.
