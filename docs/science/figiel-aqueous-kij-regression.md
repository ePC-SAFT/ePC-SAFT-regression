# Figiel Staged Aqueous Current-Catalog Recovery

Status: `READY_WAITING_PROVIDER_AQUEOUS_ACTIVE_SOLVATION_FACTOR_DERIVATIVE`.
The source and Regression contract are frozen. The existing active-Born and
aqueous-active-`k_ij` Provider seams are sufficient for stages A and C, but no
installed Provider artifact yet exposes the exact active water-solvation-factor
derivative required by stage B. No staged runtime, wheel, Validation campaign,
catalog admission, prediction, or authority transfer exists.

This document is the sole design and science owner for the staged Figiel
current-catalog recovery. It replaces the falsified assumption that one
fixed-family `164 x 11` fit should recover every printed Table 4/5 interaction.
The rejected fit remains immutable evidence; it is not proof that the
source-described staged procedure is impossible.

## Bounded claim and sequence

The smallest source-backed sequence is:

1. **Stage A — Born diameters (`5 x 5`, rank 5).** Reuse the existing first
   tracer: five active Born diameters for Li+, Na+, K+, Cl-, and Br- against
   the five SI Table S5 reported-average water-solvation targets.
2. **Stage B — water solvation factor (`21 x 1`, rank 1).** Fit one
   ion-independent `f_water` to all 21 audited Hamer--Wu NaBr MIAC rows.
3. **Stage C — aqueous interactions (`164 x 11`, rank 11).** Fit the eleven
   Table 4/5 interactions to all 164 audited Hamer--Wu rows.
4. **Confirmation.** Replay A-B-C with the preceding stage outputs fixed. Run
   at most three confirmation cycles and require the maximum scaled coordinate
   change between consecutive cycles to be at most `1e-5`.

The final candidate either reproduces the printed Table 4/5 tuple within the
frozen `0.05` maximum parameter difference and passes the in-sample observable
gates, or retains
`SOURCE_DESCRIBED_STAGED_RECOVERY_DID_NOT_REPRODUCE_PRINTED_TUPLE`. Failure
does not authorize changing rows, weights, bounds, starts, cycle limits, or
tolerances.

The eleven immutable Stage-C coordinates are:

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

The explicit fitted zero is a target. Blank Table 4/5 cells are not zeros.
Inherited water--methanol and water--ethanol interactions are excluded.

## Source statements and Regression choices

The primary paper is Figiel, Yu, and Held, *Industrial & Engineering Chemistry
Research* 64 (2025) 9406--9418, DOI `10.1021/acs.iecr.5c00475`. The retained
71,826-byte Markdown has SHA-256
`ce80533925a91bc59d8d0d8056113c40611ca26c2edf04aced76986d50bd4bae`.
Durable line locators are 277, 279, and 281:

- line 277 says Born diameters were adjusted to water solvation Gibbs energies
  and solvent-specific, ion-independent `f_k` values to NaBr MIAC data;
- line 279 says aqueous ion--water and ion--ion `k_ij` values were adjusted to
  experimental aqueous-salt MIAC literature data; and
- line 281 declares the order Born diameter, `f_k`, then `k_ij`, followed by an
  iteration whose parameter changes were small.

The paper does **not** disclose an exact MIAC objective, weights, row subset,
bounds, starts, or cycle termination. Those are Regression-owned choices below,
not reconstructed author choices.

Held et al. (2014), retained Markdown SHA-256
`b8b1e46bf870224de5de68b5989f9cb377d17445d87109a5462a94f1efaafbda`,
provides lineage evidence rather than a Figiel-specific prescription. Lines
255--268 define a squared relative-deviation objective for osmotic coefficients,
including Eq. 20, `sum(1 - phi_calc/phi_exp)^2`, and describe sequential
parameter fitting. Regression therefore uses the analogous relative MIAC
residual, while explicitly labeling that choice.

The approved Validation packet is commit
`8944d34f7002cda1bb8760e606cc1f11696f58cd`, tree
`6c8fd350dcd6bfdd7be1918f73fd33a23e2070dd`, with:

- target-ledger SHA-256
  `f405a3e48d21cd979a8dd480d5f8cb3be40754f5d6babf368b505b5f305607f0`;
- parameter-packet SHA-256
  `932e8baa90fcefbaa8c3a8730cdeadd83a4c01f0a3b109f4e4cd0319aee9312b`;
- metadata SHA-256
  `8ea06c6ca5452d01448a03f9a76cf7d0c35bb99c9abe23ccb1729d56c71d468f`;
- SI extraction SHA-256
  `85bd39f727158d5a9d6eea6828c1673f73850e783a655b09660cc9b66d84321a`;
  and
- Hamer--Wu CSV SHA-256
  `2f63e13f06a5b0f4e8bca2980b6a8d9d7fb0f839153c43e3a71952daf9796595`.

All states are aqueous molality-scale observations at `298.15 K` and
`100000 Pa`. Stage A uses exactly the five SI Table S5 reported averages
already frozen by the Born-tracer contract. Stage B uses every one of the 21
audited NaBr rows in the source packet (`0.001` through `6 mol/kg`). Stage C
uses every one of the 164 LiCl, NaCl, KCl, LiBr, NaBr, and KBr rows. Each row
appears once. The `<=5 mol/kg` subset is sensitivity evidence only. All rows are
training data; there is no pointwise uncertainty or approved held-out cutoff.

## Frozen residuals, variables, and solver

Stage A retains its existing residual, units, bounds, scales, starts, exact
Provider Jacobian, and closed result contract without reinterpretation.

For Stages B and C, positive observed molality-scale MIAC
`gamma_q_observed` defines the dimensionless equal-weight residual

```text
r_q(theta) = 1 - gamma_q_model(theta) / gamma_q_observed.
```

No observed value is treated as uncertainty. Equal weights are a predeclared
Regression choice because the packet has no rowwise uncertainty.

Stage B has one dimensionless variable `f_water = z`, scale `1`, bounds
`[1,2]`, and starts `1.2` and `1.8`. Published `1.5` is comparison-only and is
not a seed. Stage C has eleven dimensionless variables `k_j = z_j`, each scale
`1`, bounds `[-1,1]`, and starts all `0`, all `-0.5`, and all `+0.5`. The
published tuple is comparison-only and is not a seed, prior, or regularizer.

Every stage uses the existing Ceres owner with `DENSE_QR`, one thread, silent
logging, at most 500 iterations, and function, gradient, and parameter
tolerances `1e-10`. No second Ceres engine, result family, native module, or
CMake target is admitted.

## Exact Provider derivative contracts

Stage A consumes the existing model-bound active-Born value/first-total-
derivative callback.

Stage B requires one appended model-bound callback for an ordered aqueous
`(water, Na+, Br-)` model. At fixed `T`, fixed `P`, and formula-unit molality,
it accepts one finite trial `f_water` and returns
`ln(gamma_pm^m)`, exact total fixed-pressure
`d ln(gamma_pm^m)/d f_water`, reference convergence diagnostics, parameter
fingerprint, component order, and structured status. Provider retains EOS,
Born/association/electrolyte reference sequence, density closure, and CppAD
ownership. Regression derives

```text
gamma_model = exp(log_gamma_model)
dr/df_water = -(gamma_model/gamma_observed) * dlog_gamma_model/df_water.
```

Stage C consumes the existing Provider callback from correction
`06d21af0334a22bafd31d617f3c8535b53711140` (merged as
`39b39d9d7831a8da943372df019f5a9d7d388b44`). Its retained wheel SHA-256 is
`3bacac6818708091629a79ce9a7a320a07f87093f0697508060aa1aee7368cb6` and
installed-header SHA-256 is
`01568808f48c8cf0cd5fd0eb0b3d038349a319251eb8a9f9053b028ed35e5a36`.
For each ordered `(water,cation,anion)` row it returns `ln(gamma_pm^m)` and
exact total fixed-pressure derivatives with respect to
`(k_water_cation,k_water_anion,k_cation_anion)`. Regression maps those three
entries into the global eleven-column Jacobian and computes

```text
dr/dk_j = -(gamma_model/gamma_observed) * dlog_gamma_model/dk_j.
```

All other row entries are structural zero. Production numerical derivatives,
copied EOS/reference equations, independent density closure, and Equilibrium
dependencies are forbidden.

## Derivative, rank, and confirmation gates

Installed-artifact derivative checks use centered callback-value differences at
steps `1e-4` and `5e-5`. Each active column must meet

```text
abs(J_exact - J_h/2)
  <= max(1e-8, 20*abs(J_h - J_h/2), 2e-8*abs(J_exact)).
```

Finite differences are evidence only. They are not a runtime backend.

For each stage, SVD rank uses

```text
s_max * max(residual_count, parameter_count) * epsilon_binary64 * 100.
```

Required ranks are exactly `5`, `1`, and `11`. Results retain singular values,
condition number, complete/nonzero columns, active bounds, and the
least-sensitive direction. Rank failure stops; Regression does not add priors,
regularization, parameter tying, or fixed published values to manufacture
rank.

Each declared start must converge to the same stage solution within `1e-5` in
scaled infinity norm. A primary A-B-C pass is followed by at most three
confirmation A-B-C cycles. Each stage starts from the preceding cycle's value,
not a published seed. The campaign is numerically cycle-converged only when the
largest scaled change across all 17 coordinates is at most `1e-5`. Failure
after three cycles is retained without tuning.

## Retained rejected alternatives

`evidence/figiel-aqueous-kij-published-tuple-preflight.json`, SHA-256
`5bd86e332b94781112eeee0ca06765a0f084020a30af76169861bbc610d5743d`,
retains the exact isolated logarithmic-residual result. Its installed
`164 x 11` Jacobian has rank 11; all 492 declared derivative checks pass; the
three starts agree within `6.740175084729572e-11`; fitted cost is
`0.25101017331848846` versus published cost `0.4200114112320464`; and the
maximum published-parameter difference is `1.8` with five lower bounds active.

A read-only discriminator using the frozen relative residual also selected an
incompatible boundary optimum: all-164 published AARD `5.13246%` versus fitted
`4.19795%`; `<=5 mol/kg` published AARD `4.68197%` versus fitted `3.46284%`;
maximum parameter difference `1.8`. These results justify the staged source
contract. Neither establishes global impossibility or author error.

## Result and status semantics

One closed staged result retains source/artifact identities, input and
evaluated row IDs, each stage/cycle/start, all 17 fitted coordinates, published
comparisons, row residuals, objectives, gradients, rank/conditioning,
sensitivity, active bounds, Provider diagnostics, and ordered failure reasons.

Statuses remain independent:

1. `solver_converged` is reported per Ceres stage/start.
2. `numerically_converged` requires exact derivatives, required local rank,
   declared-start agreement, finite diagnostics, and confirmation-cycle
   agreement within `1e-5`.
3. `physically_valid` requires every fixed-pressure state to be Provider-usable
   with positive finite MIAC and all existing Born-tracer physical gates.
4. `workflow_valid` additionally requires exact source hashes, row membership,
   units/basis, state, salt/component order, Provider artifacts/fingerprints,
   fixed-family inputs, and complete input/evaluated/failed accounting.
5. `scientifically_valid` additionally requires every frozen in-sample
   observable gate and
   `max_j abs(k_j_fit-k_j_published) <= 0.05`.
6. `predictive_status` remains
   `NOT_ADJUDICATED_NO_APPROVED_HELD_OUT_CUTOFF`.

The `0.05` comparison is a user-approved engineering recovery criterion, not
source uncertainty or a uniqueness statement. Existing in-sample forward
checks remain pooled MIAC RMSE `<=0.17`, per-salt RMSE `<=0.35`, per-salt
maximum absolute MIAC error `<=1.25`, and first predicted MIAC for each salt
`<0.98`. A failure remains a valid scientific result but not a valid recovered
catalog candidate.

## Ownership, next gate, and negative space

Regression owns targets, residuals, Ceres execution, diagnostics, and the one
staged result. Provider owns values, exact derivatives, model records,
reference sequences, and density closure. After an exact Provider wheel/header
implements and verifies the Stage-B seam, Regression may add only the compact
164-row packet and focused staged cost/workflow owners to its existing package,
native module, result family, and target.

The next gate is therefore one Provider artifact with the approved active
water-solvation-factor callback and installed derivative proof. Runtime work
must then prove ranks `5/1/11` before retaining a staged candidate. Regression
may author the later installed-artifact Validation campaign under the standing
serialized writer rule; its own evidence cannot self-promote.

Excluded: generic registries, mutable parameter overlays, Provider catalog
writes, compatibility shims, simultaneous all-table solves, organic-solvent or
expanded-ion scope, density/osmotic targets, association/polar/reactive/MEA
scope, uncertainty, prediction, global identifiability, a second solver/result/
module/target, and any runtime dependency on Zotero, Validation source trees,
Migration, or lab code.
