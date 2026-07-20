# Pure-Saturation Propane Checkpoint A Plan

**Goal:** Extend the accepted one-engine pure-saturation workflow to one local,
source-backed propane candidate while preserving methane and ethane numerical
parity.

**Architecture:** Add one immutable propane dataset/specification record and
one closed runner choice. Reuse `PureSaturationFitResult`,
`fit_pure_saturation`, the single `_native` target, the existing lifted-volume
Ceres formulation, and only provider
`evaluate_pure_phase_parameters`. No equation, solver, result, target,
transport, or registry owner is added.

## Fixed inputs and boundaries

- Regression start: `e01a272e7d11e39db3cdc4105e0c460aceca637d`
  (tree `50b76650f1386f1d53e69760221028a51b964e39`).
- Coordination gate: Migration
  `76f9946a99ac0e000ecafd5644d4ed227f469b5e`.
- Source packet: Validation
  `7e51590757f1cb85f51df98e9fe1f88cd4255a88` (tree
  `05af9e948c786ddfcf43dba701970f1cbb6435a2`).
- Target CSV SHA-256:
  `ccd1cfa15ec44432b06cbf22316d168c61b282631c9b1e1591e497b8d48b5676`;
  packet YAML SHA-256:
  `ba31448989f565d05d63908076e836977780aa87199f208310e9b80b03f64697`;
  63-row source receipt SHA-256:
  `ed5eb703ccd3e6bb4c4cfa82ecd58c58f9da0c93ab07a204dee94d8b0ae8d081`.
- Dataset/specification:
  `glos-2004-experimental-propane-saturation-110-340-k-v1` /
  `pure-propane-saturation-lifted-volumes-v1`; exact component `propane`.
- Partition: train 150/210/270/330 K; hold out the 18 interior non-training
  rows; stress 110/340 K.
- Provider artifact: commit `4b10cb899c94687cae734980285badb224dc95e6`,
  wheel SHA-256
  `f92f79c8d6f614660e5c201b7061c9b02b5cd1a25a4ed8c8fee0b59adaabf2bf`,
  bundle `gross-2001-propane` v1, parameter fingerprint
  `sha256:9bfbc8d7789e51609945e61dbdf7a020decc8f9e31b408b0977724c7cb3e1551`.
- Start `(2.0020, 3.6184 angstrom, 208.11 K)` and molar mass
  `0.044096 kg/mol` are the provider's source-backed Gross 2001 records.
- Preserve the existing parameter bounds/scales, residuals/weights, Ceres
  controls, derivative contract, rank gates, and confirmation gate. Volume and
  reporting-pressure domains may widen only enough to contain the frozen source
  grid; they are solver domains, not accuracy tolerances.
- Source uncertainties remain immutable reporting context. They do not weight
  the fit and create no cutoff. Held-out/stress model errors remain descriptive.
- Do not change or execute the blocked binary formulation.

## Task 1: Freeze records and parity with RED tests

Modify `tests/test_records.py` and `tests/test_native_fit.py` first. Require the
exact 24 target rows and 4/18/2 partition, all pointwise source uncertainty
fields, exact packet/provider identities, the propane root export, and closed
rejection of aliases. Parameterize the existing training residual/Jacobian
directional check for propane. Before production edits, run:

```text
python -m pytest tests/test_records.py tests/test_native_fit.py -q
```

Record RED only for missing propane ownership. Keep the existing methane and
ethane expected parameters, costs, ranks, and reporting predictions unchanged.

## Task 2: Add the smallest production record and engine branch

Copy the target CSV byte-for-byte to
`src/epcsaft_regression/data/propane_saturation.csv`; retain packet commit/tree
and all three governing hashes in `src/epcsaft_regression/records.py`. Add
`PROPANE_SATURATION_FIT_V1`, a dedicated exact-header parser for the richer
Glos row schema, and optional uncertainty/vapor-density observation fields.
Do not copy the Validation YAML or 63-row audit receipt into the runtime wheel.

Extend only the closed component validation in
`src/epcsaft_regression/native/pure_saturation_contract.cpp`. The existing
`pure_saturation_fit.cpp`, provider transport, result class, Python workflow,
CMake target, and Ceres residual/Jacobian equations remain single owners.
Run record and native tests to GREEN.

## Task 3: Extend the one runner and prove candidate behavior

Add `--component propane` to `tools/run_candidate.py`, selecting
`gross-2001-propane` v1 and the immutable propane specification. Add a distinct
predictive status field if the existing result does not already expose it;
solver, numerical, physical, and predictive statuses must not be inferred from
one another. Generate one canonical propane receipt and check a second run for
byte identity before claiming runner stability.

The candidate must have finite in-bound parameters, Ceres convergence, full
11-column rank and fitted-parameter rank 3, deterministic confirmation, and
the existing numerical/physical closure gates. A scientific mismatch is
reported without inventing a model-accuracy tolerance.

## Task 4: Verify and retain one commit-bound wheel

Run the full source suite, build one wheel from the green implementation
commit, and retain it below
`../artifacts/regression-pure-saturation-v1/<implementation-short-hash>/`.
In a fresh Python 3.13 environment install only the exact provider and
regression wheels, then run the full package suite and candidate runner. Audit
wheel members, imports, native linkage, absence of source/sibling paths, one
CMake native target, three packaged data files, and unchanged binary evidence.

Update `CONTEXT.md`, `ARCHITECTURE.yaml`, `README.md`,
`docs/science/pure-saturation-regression.md`, and the current capability record
with observed results and hashes. Claim archive-byte reproducibility only if a
second clean build from the same implementation commit matches byte-for-byte.
Commit locally and stop for permanent-lab adjudication; do not push, promote,
publish, persist fitted parameters, or create a migration receipt.
