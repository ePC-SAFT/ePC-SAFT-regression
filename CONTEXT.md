# Regression Repository Context

This repository owns one strict pure-saturation regression workflow. Methane
and ethane are accepted. The propane Checkpoint A and constant methane/ethane
`kij` subjects are retained authority-neutral, non-production provenance; D-023
removes both from the active migration sequence.

`governance_doctrine_revision: 3`

Canonical local doctrine: `../ePC-SAFT-organization/GOVERNANCE.md`.

Regression owns both its package evidence and execution of assigned
installed-artifact campaigns in the sibling Validation repository. The task
remains based here; Validation remains the durable black-box evidence home and
has no resident worker. Migration serializes exact campaign subjects and
dispatches a separate reviewer only when required. This execution model does
not change scientific admission or runtime authority.

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

`parallel_regression_status: FIGIEL_BORN_TRACER_PACKAGE_CANDIDATE_BLOCKED_PUBLISHED_DIAMETER_RECOVERY`

`runtime_status: AUTHORITY_NEUTRAL_PACKAGE_IMPLEMENTED_VALIDATION_NOT_STARTED`

Migration D-023 supersedes D-022's active order without rewriting its evidence.
Canonical Migration checkpoint
`3a4ef0a0c6b98c43405d3cafc1ac4f5f87afa68d`, tree
`9307c3f79581b6e0479d4ac2468932b2a68e5f5b`, preserves the published Figiel
design as a useful parallel, deferred Provider-dependent track while HELD2
retains Provider implementation priority. The frozen neutral-hydrocarbon
equations and numerical contract remain in
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

The source-faithful Figiel, Yu, and Held 2025 current-catalog ion Born-diameter
tracer is a complete, published design at Regression
`8191dcc9fc038caac1f52cd22303c600e2b61241`, tree
`656cc22409e1e1632f536184e2f719640178748f`, retained in parallel behind the
HELD2 critical path. Provider design
`da9660481f08bb5557cc03da528edec15cc15e1f`, tree
`e34575ae646c40f3fb63b5994c957e30bb035f69`, and Validation source packet
`8944d34f7002cda1bb8760e606cc1f11696f58cd`, tree
`6c8fd350dcd6bfdd7be1918f73fd33a23e2070dd`, are permanent-lab approved.
The complete source-bound Ceres contract is singularly owned by
`docs/science/figiel-born-diameter-tracer.md`. It fixes exactly five Table S5
reported-average targets, five active Born diameters, the diagonal exact
Jacobian, source-resolution comparisons, confirmation, and statuses.

Migration D-027 bound Provider implementation
`907b077ec6f841a8a028fc759df14f899c79339c`, wheel SHA-256
`c327b9a176e54bfc79b625cca7f0c87f2a62fc7d87059826e40c9d70e214f0cd`,
and installed-header SHA-256
`610cc480f05c3e17e431d26fd1b2c8628eec3e2adb412102a284d4d5d6eb8171`.
The package implementation consumes that callback through the existing native
module and target. All three starts converge, every scaled residual is below
`1e-8`, the `5 x 5` exact Jacobian has rank 5, and no bound is active. The
source observables round-trip, but fitted-minus-published diameters range from
`-0.003001472493401991` to `+0.007461646407642686 angstrom`, exceeding the
frozen `0.0005 angstrom` comparison. The candidate is therefore scientifically
non-admitted without changing any gate. Validation writing has not started.

`electrolyte_born_parameters` remains `NOT_READY`: this package candidate must
not be interpreted as downstream readiness, predictive evidence, or Provider
catalog authority.

The eventual direction is a staged, family-by-family recovery campaign for
Figiel Tables 2–5. That later intent does not broaden the first Born-diameter
tracer into a simultaneous all-table fit or generic registry. Validation's
future exact ledger alone will classify each value as fitted, inherited,
fixed, or blank; Regression makes none of those assignments in advance.
