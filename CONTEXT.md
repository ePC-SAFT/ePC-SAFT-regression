# Regression Repository Context

This repository owns ePC-SAFT target contracts, Ceres parameter fitting,
diagnostics, and results. Methane and ethane pure-saturation workflows are the
accepted production capabilities. Other completed campaigns remain
authority-neutral evidence unless separately admitted.

The active Python package is version `0.2.0.dev0` and consumes exactly
`epcsaft==0.2.0.dev0`. Public workflows construct EOS `Parameters` and
`Mixture` objects; source-bound trial models use `Parameters.from_records`.

`provider_frontend_status: PARAMETERS_MIXTURE_0_2`

The active installed EOS subject is commit
`24ab1bdfb3a1558eb87815e9c3b5d97fbe6a025d`, tree
`145cf8b6126ad403aa8ee51b1b4623db3da8aa2f`, wheel SHA-256
`66b7ea8fb29e0a268b555cbdf401c3502517c088669a4157e8f64ab985b59ce9`, and
installed public-header SHA-256
`db37805c1abd9b1a355f41a89154d14756477c7293af940dc298a6b038aee45d`.
Its installed static-library SHA-256 is
`5e5c3e8311c011365471eeff29b2d9ff8486d149745a4bea7cd48bd6e784f9d3`.
The frontend migration changes public construction and generic capability
names only; it does not change accepted authority, target contracts, or
numerical acceptance criteria. Retained pre-rename Provider 0.1 receipt
identities remain immutable historical evidence and are not active replay
inputs.
The previous frontend subject's final EOS 0.2 JUnit receipt artifact is
`provider-python-frontend-v0.2/14fa374/provider-tests.xml`,
SHA-256
`56d1d3e9fcf47bb700df4223fa2d9a20444dc97aa22b5c39525d446a296ba3cf`, with
276 passed and 1 skipped tests. It remains historical evidence for that exact
subject and is not rebound to the active ordinary-sigma artifact. Portable
artifacts are resolved relative to an executor-supplied immutable artifact
root; workstation paths are execution metadata only.

`governance_doctrine_revision: 5`

Canonical doctrine:
`https://github.com/ePC-SAFT/ePC-SAFT-management/blob/main/GOVERNANCE.md`.
During the local checkout transition, resolve doctrine and the access manifest
from `../ePC-SAFT-organization/GOVERNANCE.md` and
`../ePC-SAFT-organization/manifests/management.json`.

EOS is the canonical active repository and role name. Retain `Provider` only
inside historical receipt prose, artifact paths, schema fields such as
`provider_frontend_status`, and public ABI/API symbols whose spelling is part
of an accepted contract.

Regression owns both its package evidence and execution of assigned
installed-artifact campaigns in the sibling Validation repository. The task
remains based here; Validation remains the durable black-box evidence home and
has no resident worker. Management serializes exact campaign subjects and
dispatches a separate reviewer only when required. This execution model does
not change scientific admission or runtime authority.

Literature work follows
`docs/science/literature-reproduction-contract.md`. Regression retains compact
routine numerical-tolerance sentinels for its generic public workflows and
supported active-coordinate variants. Validation owns full literature data
grids/objectives, broad benchmark matrices, profile and uncertainty campaigns,
induced/cross-association comparisons, transfer evidence, and durable
installed-artifact receipts.

Accepted historical receipts `promotion-0020-regression-methane-saturation-v1`
and `promotion-0023-regression-pure-saturation-ethane-v1` make this repository
the production owner of the exact reproducible methane and ethane workflows.
Validation receipt `validation-0022-regression-pure-saturation-ethane-v1`
passed, and state receipt `state-0025-regression-ethane-publication` verifies
ethane publication. Neither workflow gives its fitted parameters predictive,
uncertainty, scientific, or EOS-catalog admission.

Both components use the installed EOS SDK, one Ceres engine, exact
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

`parallel_regression_status: FIGIEL_BORN_TRACER_PACKAGE_CANDIDATE_PASSED_LOCAL_GATES_VALIDATION_NOT_STARTED`

`runtime_status: AUTHORITY_NEUTRAL_PACKAGE_IMPLEMENTED_VALIDATION_NOT_STARTED`

`next_figiel_family_status: CONDITIONAL_STAGE_C_NUMERICALLY_CONFIRMED_PRINTED_TUPLE_NOT_RECOVERED`

`general_parameter_regression_status: TEN_PARAMETER_FAMILIES_FIT_READY_AUTHORITY_NEUTRAL`

`multi_parameter_core_status: IMPLEMENTED_ACCEPTED_PURE_WORKFLOWS_CUT_OVER_LOCAL`

The canonical general design and native core now implement one ordered
`[N fitted | Q lifted]` Ceres layout, complete ordered fitted-parameter start
vectors, observation-owned lifted starts, an
explicit evaluator-slot sharing map, exact full and nuisance-projected
Jacobians, and residual-only versus exact-Jacobian request propagation. A compact
closed-form native test proves `N = 2`, `Q = 1`, full rank 3, projected rank 2,
fitted-vector confirmation, explicit slot sharing, structural rejection,
incomplete-buffer rejection, and missing-column rank failure without finite
differences.

The installed EOS's joint pure callback is the first real `N > 1` seam.
Accepted methane and ethane training now use the shared contiguous core with
`N = 3`, `Q = 8`, `R = 16`, full rank 11, and nuisance-projected parameter
rank 3. Their existing public presentation wrapper and common-pressure
reporting solve remain, but the former duplicate pure training Ceres loop is
gone. Frozen methane/ethane fitted values, costs, reporting predictions, and
statuses remain within their existing numerical contracts. This does not
imply that arbitrary simultaneous mixture parameters are executable; each
such block still needs one installed exact multi-active evaluator contract.
The public `fit_parameters` transport admits the same closed ordered
`(m,sigma,epsilon/k)` pure block and returns three ordered parameter
diagnostics. It still rejects every other arbitrary `N > 1` block; the pure
adapter is not a generic cross-package callback.
The installed pure callback returns value, gradient, and Hessian together;
residual-only Ceres calls avoid copying a Regression Jacobian but cannot avoid
the EOS derivative computation. That artifact limitation is not reported
as a sensitivity-free value-only evaluation.

The campaign-specific EOS subjects and artifact hashes recorded below are
retained historical evidence for their named reference campaigns. The active
installed frontend binding for new replay is the `provider_frontend` subject
above; historical evidence is not rewritten as if it had used that artifact.

The installed EOS advertises one ordinary-`sigma_ij^3*kappa`
fixed-topology association descriptor and variable-width callback. Regression
selects typed ordinary-pure and component/site-pair slots, keeps `N` fitted
coordinates distinct from `S` EOS slots, and accumulates exact Hessian columns
for explicit sharing. Vapor-pressure, fixed-pressure-density, and combined
saturation rows use this one path. A 2B five-parameter case and a non-2B
multi-pair case prove mechanics only; Regression packages no MAPA/DEEA dataset,
accepted topology, or accepted tuple. Status is
`FIXED_TOPOLOGY_ASSOCIATION_MECHANICS_EXECUTABLE_NO_PARAMETER_ACCEPTANCE`.

The user-approved 2026-07-27 destination is general parameter-family
regression, defined in `docs/science/general-parameter-regression.md`. The
first ten families are implemented. Caller-supplied, source-bound
fixed-composition VLE rows may independently fit one shared `k_ij` or `l_ij`
for any neutral, nonassociating binary model whose installed EOS
descriptor advertises the corresponding exact `(n1,n2,V,pair_parameter)`
Hessian callback. The runtime uses the existing native module
and target, one Ceres engine/result, exact derivatives, caller-declared
bounds/scales/starts/volume policy, full and projected rank diagnostics, and
isolated training/reporting partitions. Unknown EOS descriptor members
are reported unsupported and cannot be requested dynamically.

The ninth family is model-level
`ionic_region_relative_permittivity`. On an installed EOS advertising
the exact SSM+DS solvation-Gibbs value/first derivative, callers may fit one
coordinate from source-bound single-ion solvation targets through the same
direct-observable Ceres path. Validation subject `81e9a3f` independently
fits all five admitted ion models against the Table S5 reported averages
while holding the paper's Born diameters fixed. Each `1 x 1` fit is full-rank,
non-bound, confirmed from starts 4 and 12, and returns
`7.99999999787946` through `7.999999997884154`. Figiel reports `8` as fixed
model input, so this establishes fit-ready mechanics rather than recovery of
a paper-fitted parameter or a shared multi-model parameter.

The tenth family is component `relative_permittivity` on an installed
single-ion solvation callback that identifies the active solvent separately
from the observed ion. EOS subject `7cadaad` supplies the exact total
first derivative through bulk dielectric, Debye--Huckel, Born,
reference-state, and fixed-pressure density paths. Regression subject
`177890a` reuses the same direct-observable Ceres owner. Corrected Validation
subject `e4cb7af` fits water independently against all five Table S5
reported-average solvation targets with the paper's Born diameters fixed.
Every `1 x 1` problem is full-rank, non-bound, and confirmed from starts 50
and 110, returning `78.0899937514462` through `78.08999375166104`. The fixed
input `78.09` and its `0.005` reporting half-increment are descriptive
comparison evidence, not uncertainty or scientific acceptance.

Polar parameter families are outside the Regression roadmap by user decision.
They remain representable EOS inputs where applicable, but no Regression
derivative, residual, fit, or validation work is planned for them. This does
not relabel them unsupported EOS physics.

The coordinated MEA path is owned by
`docs/science/mea-coupled-regression-master-plan.md`. The first coupled tracer
is one admitted CO2 partial-pressure observation plus one eligible liquid
speciation equality or aggregate with at most two fitted coordinates. Loading
is a state input; heat is not admitted because no source-complete heat contract
exists. The first tracer uses the user-approved low-pressure
`p_CO2 = f_CO2_liquid` convention from a certified homogeneous reacting
liquid; it needs no vapor composition or bubble solve. EOS exact active
parameter partials and Equilibrium conditioned active-parameter amount/volume
sensitivities are implemented. The generic process-local positive-observation
adapter is implemented over the existing Ceres core with exact value/Jacobian
request propagation and fail-closed metadata and certificate checks. The
retained installed Hilliard/Böttinger state is `FEASIBLE_ONLY`: its solver,
numerical, physical, and EOS-domain gates pass, but an independent
reduced-Hessian calculation confirms negative curvature and the local-minimum
gate fails. A certified stable installed state and frozen rank-sufficient
`N <= 2` coordinate block remain prerequisites. Reactive bubble issue 37 is
closed as not planned and is not a tracer gate.
The application
candidate is a staged three-coordinate ion block, not the stale historical
12-parameter block. The frozen pressure/speciation split has at most 297
eligible training and 267 eligible reserved residuals.

Acquiring the missing Gross--Sadowski/Cameretti pure-association series or a
standalone cross-association recovery dataset remains useful optional evidence,
but it is not the MEA critical path.

Caller-supplied pure-saturation pressure and liquid-density rows may
independently fit one component `segment_count`, `segment_diameter`, or
`dispersion_energy_over_k` coordinate for a neutral, nonassociating,
constant-diameter pure model. The same general Ceres owner uses two lifted log
volumes per row and four scaled residuals: liquid and vapor pressure closure,
chemical-potential equality, and liquid-density reproduction. Its exact
Jacobian consumes the EOS `(n,V,active_parameter)` Hessian, reports full
and nuisance-projected rank, checks mechanical stability and phase-volume
ordering, and reruns every declared start. This independent-family surface
does not replace the accepted joint `(m,sigma,epsilon/k)` workflow.

EOS commit `1e571ab0a84603a51ed6994b14286f683fb12b88` retains those
compatibility and parity corrections and appends the exact active-`l_ij`
capability beside active `k_ij`. The retained EOS wheel has SHA-256
`6536edc63adaa13c5c6c67c185d82c9ae232048e99dc3dc3be502708eea4410f`;
its installed public header has SHA-256
`b667379c2d7106d012c6b57f96b6f32dd23ef305fe6f15a87c22ab20029008f8`.
Regression commit `da7a44ce093201022aec2f3514d4e4fd9d8d2929` extends the same
bounded public surface and Ceres owner to the distinct `l_ij` identity.
The retained Regression wheel has SHA-256
`3a59d2233fec51f949a7784937b54f7f66beae2476fa8e33976672a480b67137`.

EOS commit `86983ff` (tree
`cecdfddde1c29f6d33fc1106d0662e0376981e4a`) appends the scalar pure
capabilities without another EOS tape. Its retained wheel SHA-256 is
`9df62965f55876104585504f0f9170fefffbbaa9e94b23f90ac1da582cf5cb4c`
and installed header SHA-256 is
`feadc42414188cd9fa50708dff5e7815f9eec85bf868f080e65e07a4c70b461d`.
The four accepted methane training rows give full lifted rank 9 and projected
parameter rank 1 for each independent fit. The fitted values are
`m = 1.0001569260577763`, `sigma = 3.7063548743836034 angstrom`, and
`epsilon/k = 150.00325287725062 K`, all non-bound with confirmation-start
agreement. These are local in-sample implementation anchors, not fitted
catalog values or predictive acceptance.

EOS commit `2d1816cf376294156684fee85611a93fc41d0970` (tree
`0a36f875fe75cc90cece5144373ab031c8724ccc`) advertises the existing
active-Born and aqueous-solvation-factor callbacks through the same typed
capability negotiation, including model-topology fingerprints. Its exact wheel
SHA-256 is
`bc40a9ac2c217afe163cf0bed36159e8f3387d6a1bd09a24665b4a952d84a8d3`
and installed-header SHA-256 is
`568381308595713a60ee0b24b1ca357956b51b7cbbf759cc13404108eeac8780`.
The shared Regression surface now fits one Born diameter from source-bound
solvation-Gibbs observations or one solvation factor from source-bound mean-
ionic-activity observations. The Figiel reference evidence consists of five
independent `1x1`, rank-1 Born problems and the existing `21x1`, rank-1 NaBr
water-factor problem. All starts converge to the prior paper-specific
reference values with no active bounds. This is fit readiness only for the
advertised fixed-state reference paths, not general electrolyte, catalog, or
predictive authority.

EOS commit `e063f32f6e2b975af11ae7e65c1a12cd3076941c` exposes the
model-bound Figiel dielectric ion-suppression coefficient through an exact
relative-permittivity value/first-derivative capability. Its retained wheel
SHA-256 is
`d090fa92d86ea91ab40df24aa96fc18c1821eaf636e9a03e01346a6bef8e31bd`;
the installed public-header SHA-256 is
`5dccb80d0d64a24bb8070bd5a0d7f27c3b55173253b1637f9a8a7e03e4e6f82c`.
Regression commit `77a94c215d3dadb52a4f0ed90968631949b612f4`
fits one coefficient from caller-supplied salt-free-normalized relative-
permittivity rows through the existing engine/result/native target. The
36-row Figiel Figures 2-3 water/methanol reference campaign is rank 1,
non-bound, and converges from both declared starts to
`7.067350349980952`, compared descriptively with the paper's `7.01`.
These digitized training rows establish fit-ready execution only.

EOS commit `621adbf9f60c75a1f108256c51eb1ae14161d6b9`
fail-closed advertises exact first derivatives of infinite-dilution
ion-solvation Gibbs energy with respect to the three pair coordinates for the
source-bound methanol/K+/Br- and ethanol/Na+/Cl- selections. Its wheel SHA-256
is `8c942e9177bec1502345b5193642460166d27fcee840663a18d19d7ea3d823e2`;
the installed header SHA-256 remains
`16c1996e592808ed2764de7b17bead6ed06789182fe8c0e649e049e83fe9a348`.
The shared Regression path independently fits one advertised pair from one or
more caller-supplied solvation-Gibbs endpoints. The four retained Figiel
near-pure endpoints are constructed, rank-1 implementation evidence. Their
differences from the rounded Table 5 values are descriptive and do not claim
exact recovery of the paper's unavailable pure-organic fitting targets.
Regression subject `e9ce5a6c4441cd8ecc6cd18be2bb4d8614e741a1`
produced wheel SHA-256
`efc3b421ad94f4448ffce630528911363cc5be7565de2a5fca0e7787a6e68f2d`.
Validation subject `9857a1ac058038f8b7664cfe41ab06729f661715`
retains the byte-stable installed-artifact record with SHA-256
`e34f0033de804afcc3390002dee2e153a1dd80a42b0fe01080d965257a1ba0f9`.

All 17 audited May rows reproduce the frozen `68 x 35` formulation with
convergence, full rank 35, projected parameter rank 1, a non-bound
`k_ij = -0.00843032298906253`, and declared-start confirmation. The same
rows independently fit `l_ij = -0.002774426668544412` with the same ranks,
no active bound, and two confirmation starts. This remains in-sample,
authority-neutral evidence. It does not reinterpret the retained negative
pressure-closure evidence, establish prediction, or admit fitted values to a
EOS catalog.

Published methane, ethane, Figiel, and later campaigns remain positive or
negative evidence rather than hard-coded mixture limits. The finite typed
schema is not an arbitrary residual plugin or mutable registry. All parameter
families other than these exact neutral-binary fixed-composition, scalar pure-
saturation, advertised ion-solvation, solvation-Gibbs, and mean-ionic-activity
contracts retain their prior readiness states.

Historical Migration decision D-023 superseded D-022's then-active order
without rewriting its evidence. The retained Migration checkpoint
`3a4ef0a0c6b98c43405d3cafc1ac4f5f87afa68d`, tree
`9307c3f79581b6e0479d4ac2468932b2a68e5f5b`, preserves the published Figiel
design as a useful parallel, deferred EOS-dependent track while HELD2
retained EOS implementation priority. The frozen neutral-hydrocarbon
equations and numerical contract remain in
`docs/science/neutral-hydrocarbon-next-slice.md`; its rejected replay harness
and packaged planning framework remain absent. The new source-bound general
runtime does not promote or reinterpret that old evidence. Validation's direct-
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
resolution limit, not a single-ULP pressure plateau or an EOS-defect
finding. Receipt `evidence/propane-candidate-fit-receipt.json` has subject
SHA-256 `c059e5381a3fbf0cf00a43aaa0cc28b67c074e325f1a33066b6937c16cc0a761`
and file SHA-256
`daaa93a2f5d0e0dbf81fc83562bb809c687a36e66578071c68ea66b1988ccaec`.
No source uncertainty is treated as an acceptance cutoff, and predictive
status remains `NOT_ADJUDICATED_NO_APPROVED_HELD_OUT_CUTOFF`.

Historical review verdict `CORRECTION_DESIGN_JUSTIFIED` produced a design and
evidence protocol for
`abs(P_phase - P_report) <= atol_resolution + 1e-8*abs(P_observed)`.
The relative term is unchanged, observed pressure remains only a magnitude
scale, and no value of `atol_resolution` is selected. The derivation and replay
protocol is frozen in the two existing science owners at documentation commit
`f8eee367fb54295c90cd9d5e7d8e8a73e4b8a1ae`, tree
`e60cb9262eaef886479beef842549650d61cdf10`. Under D-023 it is provenance only:
no numeric selection, runtime correction, re-adjudication, or active blocker
follows from it.

Future Regression checkpoints use exact-subject independent subagent review;
permanent-lab approval is no longer an execution or authority gate.

For binary `kij`, EOS derivatives, the `68 x 35` Jacobian, rank 35,
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
HELD2 critical path. Historical Provider design
`da9660481f08bb5557cc03da528edec15cc15e1f`, tree
`e34575ae646c40f3fb63b5994c957e30bb035f69`, and Validation source packet
`8944d34f7002cda1bb8760e606cc1f11696f58cd`, tree
`6c8fd350dcd6bfdd7be1918f73fd33a23e2070dd`, are immutable reviewed inputs.
The complete source-bound Ceres contract is singularly owned by
`docs/science/figiel-born-diameter-tracer.md`. It fixes exactly five Table S5
reported-average targets, five active Born diameters, the diagonal exact
Jacobian, source-resolution comparisons, confirmation, and statuses.

Historical Migration decision D-027 bound EOS implementation
`907b077ec6f841a8a028fc759df14f899c79339c`, wheel SHA-256
`c327b9a176e54bfc79b625cca7f0c87f2a62fc7d87059826e40c9d70e214f0cd`,
and installed-header SHA-256
`610cc480f05c3e17e431d26fd1b2c8628eec3e2adb412102a284d4d5d6eb8171`.
The package implementation consumes that callback through the existing native
module and target. All three starts converge, every scaled residual is below
`1e-8`, the `5 x 5` exact Jacobian has rank 5, and no bound is active. The
source observables round-trip. Fitted-minus-published diameters range from
`-0.003001472493401991` to `+0.007461646407642686 angstrom`. On 2026-07-21 the
user corrected the contract: the published Table 3 diameters were not residual
targets, and their `0.0005 angstrom` reporting half-increment is descriptive
source-resolution context rather than a scientific veto. The candidate now
passes its local solver, numerical, workflow, and source-observable gates;
Validation writing has not started.
The blocked-gate package evidence at Regression `2df4a305` remains immutable
Git provenance. Corrected implementation commit
`4eb122d4e40e08fb5e9cea94bfb2fe25ccfcab79`, tree
`85d1f980ce986e1ab48ed417704fd74c95c68b66`, retains the commit-bound wheel
with SHA-256
`2210b172d60141c27dec4cc7a92c9cce812bbeda75018d81f5fcef7b15e287e5`.
Canonical package evidence is
`evidence/figiel-born-diameter-candidate.json`, file SHA-256
`99d46eafbdae3428f690543364096fb414b818db201d1c35b0c0da8b03ae91d5`,
subject SHA-256
`55ea2cd69af62c45b26179cfab6939760de23058b5a7e8c880a79f67faa417ed`.
The commit-bound wheel passes the 49 non-Born tests against the retained
accepted EOS runtime. Against assigned EOS subject `907b077`, all six Born
tests pass but three legacy numerical anchors do not; that compatibility result
is reported without moving their accepted expected values.

`electrolyte_born_parameters` remains `NOT_READY`: this package candidate must
not be interpreted as downstream readiness, predictive evidence, or EOS
catalog authority.

The active source-bound family is the Figiel aqueous recovery in
`docs/science/figiel-aqueous-kij-regression.md`. The first standalone
water-solvation-factor checkpoint now fits all 21 audited NaBr MIAC rows as a
`21x1`, rank-1 Ceres problem using the predeclared equal-weight residual
`1-gamma_model/gamma_observed`. Both declared starts converge, all rows are
evaluated, and the fitted value is `1.5590515389548207`. The paper's
one-decimal `1.5` is a sufficiently close descriptive comparison under the
user's engineering judgment; it is not a second numeric acceptance gate.
Every row is training evidence and predictive status remains
`NOT_ADJUDICATED_NO_APPROVED_HELD_OUT_CUTOFF`.

Historical Provider artifact `fa766e5` supplies the exact installed NaBr batch
value and water-solvation-factor derivative used by the retained Stage-B
result. Bounded historical Provider artifact `8af6e54` supplies the retained
Stage-C exact installed batch derivatives and cancellation contract. Active
replays bind the `provider_frontend` subject above. Stage C fits each of the
eleven Table 4/5 interactions conditionally
against all affected Hamer--Wu rows while holding the other ten at the printed
tuple. Both starts converge and agree within `1.3809200687386891e-6`; the
assembled Jacobian remains full rank. Nine coordinates miss the user-approved
`0.05` comparison, with worst Li+--Br- difference `0.26079436319380545`;
Li+--Cl- reaches the upper bound.
Status is
`SOURCE_DESCRIBED_STAGED_RECOVERY_DID_NOT_REPRODUCE_PRINTED_TUPLE`. No rows,
weights, bounds, starts, or tolerances were tuned to the answer.

The authority-neutral general API now consumes the same installed batch
derivative as three closed scalar aqueous capabilities: water--cation,
water--anion, and cation--anion. A caller can fit one unordered `k_ij` from
source-bound aqueous MIAC rows while supplying the other two interactions
explicitly. This expands the observation domain of `k_ij`; it does not
reinterpret the retained Stage-C comparison, infer catalog defaults, or admit
a simultaneous electrolyte parameter overlay.

The Figiel Tables 2–5 recovery remains a staged reference campaign. It does
not define the component domain of a completed parameter family. Validation's
exact ledger alone classifies each value as fitted, inherited, fixed, or
blank; Regression makes none of those assignments in advance.
