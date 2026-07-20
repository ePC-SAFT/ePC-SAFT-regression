# Regression Repository Context

This repository owns one strict pure-saturation regression workflow. Methane
and ethane are accepted. The propane Checkpoint A and constant methane/ethane
`kij` subjects are retained authority-neutral, non-production provenance; D-023
removes both from the active migration sequence.

`governance_doctrine_revision: 2`

Canonical local doctrine: `../ePC-SAFT-organization/GOVERNANCE.md`.

Accepted migration receipts `promotion-0020-regression-methane-saturation-v1`
and `promotion-0023-regression-pure-saturation-ethane-v1` make this repository
the production owner of the exact reproducible methane and ethane workflows.
Validation receipt `validation-0022-regression-pure-saturation-ethane-v1`
passed, and state receipt `state-0025-regression-ethane-publication` verifies
ethane publication. Neither workflow gives its fitted parameters predictive,
uncertainty, scientific, or provider-catalog admission.

Both components use the installed provider capsule, one Ceres engine, exact
residual Jacobians, and immutable source and specification records. The ethane
campaign fits 140, 180, 220, and 260 K; holds out 120, 160, 200, and 240 K; and
reports 100 and 280 K as domain-stress rows. The failed 100 K reporting closure
cannot veto or establish fit acceptance.

Held-out and stress errors remain descriptive because no admission cutoff was
approved. Validation therefore records
`NOT_ADJUDICATED_NO_APPROVED_HELD_OUT_CUTOFF` and retains the ethane 100 K row
as `ETHANE_100_K_EXCLUDED_FAILURE`. Final validation evidence is at
`ePC-SAFT/ePC-SAFT-validation@5a678beff38717478fd333c65e77f005cc2f6b15`,
`results/consumer-slice-2-validation-record.json`, SHA-256
`239c84788f75f8c66240c83e4f5874f112e1197dafad6273e1c8ec4efe994d24`.
The reporting-block directional Jacobian remains an explicit evidence limit,
with no added runtime test seam.

`runtime_source_of_truth: accepted-methane-and-ethane-workflows`

`deferred_propane_evidence: pure-propane-saturation-parameter-candidate-v1`

`deferred_propane_status: BLOCKED_CHECKPOINT_A_120_K_PRESSURE_CLOSURE`

`active_regression_status: WAITING_FOR_FIGIEL_SOURCE_AND_PROVIDER_DESIGN`

Migration D-023 at control-plane commit
`fab44bd6d795afd923ffe85cc8393e7bfe576c65` supersedes D-022's active order
without rewriting its evidence. The frozen neutral-hydrocarbon equations and
numerical contract remain in
`docs/science/neutral-hydrocarbon-next-slice.md`; no packaged planning
framework, binary executable owner, or runtime API exists. Validation's direct-
experimental propane packet remains accepted source evidence at
`7e51590757f1cb85f51df98e9fe1f88cd4255a88`, but it no longer authorizes active
Checkpoint A work.

The exact blocked propane evidence subject is Regression
`33228253383ab8755384170b3caac7b38733e157`; its receipt subject is
`c059e5381a3fbf0cf00a43aaa0cc28b67c074e325f1a33066b6937c16cc0a761`.
Checkpoint A implementation commit
`aab87ebd4a40cb29f21486e06687c10eb1e44624` (tree
`462cda550cf5461e4b13cc7e40630708dd336bd7`) retains one wheel with SHA-256
`32b815fc00241516f13574594af3cf631f8fe30629bbec45ac3828d9357d705f`.
The exact installed-artifact result is solver-converged and numerically
confirmed with full rank 11 and fitted-parameter rank 3, but is not physically
valid. Held-out row `glos2004-propane-sat-120-k` has liquid-pressure closure
`1.0540036887718429e-7 Pa`, or `3.293761527412009e-8` scaled by the observed
`3.2 Pa`, above the frozen `1e-8` gate. The exact installed-callback diagnostic
reproduced that pressure and residual bit-for-bit. The required continuous
volume correction is `0.42624401466815054` local liquid-volume ULP, while both
adjacent representable volumes change pressure; this is a measured binary64
resolution limit, not a single-ULP pressure plateau or a provider-defect
finding. Receipt `evidence/propane-candidate-fit-receipt.json` has subject
SHA-256 `c059e5381a3fbf0cf00a43aaa0cc28b67c074e325f1a33066b6937c16cc0a761`
and file SHA-256
`daaa93a2f5d0e0dbf81fc83562bb809c687a36e66578071c68ea66b1988ccaec`.
No source uncertainty is treated as an acceptance cutoff, and predictive
status remains `NOT_ADJUDICATED_NO_APPROVED_HELD_OUT_CUTOFF`.

Permanent-lab verdict `CORRECTION_DESIGN_JUSTIFIED` produced a design and
evidence protocol for
`abs(P_phase - P_report) <= atol_resolution + 1e-8*abs(P_observed)`.
The relative term is unchanged, observed pressure remains only a magnitude
scale, and no value of `atol_resolution` is selected. The derivation and replay
protocol is frozen in the two existing science owners at documentation commit
`f8eee367fb54295c90cd9d5e7d8e8a73e4b8a1ae`, tree
`e60cb9262eaef886479beef842549650d61cdf10`. Under D-023 it is provenance only:
no numeric selection, runtime correction, re-adjudication, or active blocker
follows from it.

For binary `kij`, provider derivatives, the `68 x 35` Jacobian, rank 35,
projected parameter rank 1, a finite non-bound result, and both perturbed starts
passed. The three declared converged equal-weight starts did not meet the
`1e-8` pressure-closure gate: row `may2015-ch4-c2h6-002` reached
`2.2691483584976503e-2`. This does not establish global infeasibility or
model/data incompatibility. No chemical-potential equality cutoff is approved,
so the pressure-resolution protocol cannot by itself make the binary
checkpoint ready. Compact evidence remains at
`evidence/binary-kij-preflight.json`; the rejected executable alternative is
immutable at Regression `47a2a8d9579a01347df7ccaa977337ad7d4047af`, receipt
subject `d51c9f0713b6a7355be719b6843e4459f41d46d16d973668d694715d36b63676`,
compact-record head `e01a272e7d11e39db3cdc4105e0c460aceca637d`, and Migration
`4fd0e6eff5fd17399573132430d13c7e34626b98`. All 17 rows are training data.
The formulation is deferred and non-production; it is not an active blocker.

The next active objective is a source-faithful Figiel, Yu, and Held 2025
current-catalog ion Born-diameter tracer. Regression may design its Ceres
formulation only after Validation returns a stable reviewed target packet from
the main paper and official supporting information and Provider returns a
permanent-lab-approved, model-bound active-Born exact-derivative seam. The
later bounded design must then freeze the source/specification contract,
bounds and scaling, exact residual Jacobian, rank and conditioning gates,
perturbed-start confirmation, separate solver/numerical/physical statuses, and
installed-artifact validation plan. Until both upstream designs exist,
`electrolyte_born_parameters` is `NOT_READY`; no target equation or provider
API is inferred here.
