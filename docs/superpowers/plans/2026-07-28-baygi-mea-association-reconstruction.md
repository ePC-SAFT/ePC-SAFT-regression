# Baygi MEA Association Reconstruction Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fit the five pure-component PC-SAFT parameters of 2B
monoethanolamine to the frozen Baygi--Pahlavanzadeh saturation correlations
with a declared smooth-L1 approximation to the paper's Eq. 8 observable
objective and exact Provider derivatives.

**Architecture:** Append one model-bound Provider callback over
`(n,V,m,sigma,epsilon/k,epsilon_assoc/k,kappa_assoc)` using Baygi Eq. 4's
explicit `d_ij^3` association basis. Generalize the existing Regression
pure-saturation owner from fixed four-row/three-parameter arrays, eliminate
each MEA row's two equilibrium volumes with exact implicit sensitivities, and
retain its one Ceres engine, result family, Python workflow, native module,
and native target.

**Tech Stack:** C++17, CppAD, Ceres, Eigen, Python 3.13, CMake, pytest.

## Global Constraints

- Baygi 2B is a 15-row correlation reconstruction, not an exact author-run
  replay or experimental-data fit.
- The local result is a Baygi Eq. 4-inspired refit candidate. Baygi's printed
  tuple does not reproduce its printed AADs under the implemented convention,
  so no source-faithful/reference-validation claim is allowed.
- Preserve all methane, ethane, and propane behavior.
- Consume exact Provider value/gradient/Hessian output; no copied EOS or
  production finite differences.
- Do not add electrolyte, reaction, binary, cross-association, uncertainty,
  catalog persistence, a second fitter, or a generic parameter registry.
- Keep solver, numerical, physical/workflow, scientific-comparison, and
  predictive statuses separate.

---

### Task 1: Provider joint associating pure-phase callback

**Files:**
- Modify: `../ePC-SAFT/.worktrees/regression-integration/src/epcsaft/include/epcsaft/native_sdk_v1.h`
- Modify: `../ePC-SAFT/.worktrees/regression-integration/src/epcsaft/_native/derivatives.hpp`
- Modify: `../ePC-SAFT/.worktrees/regression-integration/src/epcsaft/_native/phase_block.hpp`
- Modify: `../ePC-SAFT/.worktrees/regression-integration/src/epcsaft/_native/eos.hpp`
- Modify: `../ePC-SAFT/.worktrees/regression-integration/src/epcsaft/_native/eos.cpp`
- Modify: `../ePC-SAFT/.worktrees/regression-integration/src/epcsaft/_native/native_sdk.cpp`
- Test: `../ePC-SAFT/.worktrees/regression-integration/tests/test_native_sdk.py`

**Interfaces:**
- Consumes: one-component neutral 2B `ParameterBundle.from_records` model.
- Produces:
  `evaluate_diameter_basis_associating_pure_phase_parameters(model,T,n,V,m,sigma,epsilon_k,epsilon_assoc_k,kappa_assoc,result)`
  with seven-coordinate gradient and Hessian, pressure, `mu/RT`, parameter
  fingerprint, and topology fingerprint.

- [x] Add a ctypes ABI test that requires the appended callback and verifies
  its coordinate order, finite `7 x 7` derivatives, source/topology
  fingerprints, and unsupported-model behavior.
- [x] Run the focused test and retain the expected missing-tail failure.
- [x] Extend the existing CppAD pure-phase tape so all five parameter
  coordinates are independent together; expose the result through one
  append-only SDK callback without adding target or fitting policy.
- [x] Run Provider native SDK tests plus directional value/gradient/Hessian
  checks for representative MEA liquid and vapor states.
- [x] Commit the green Provider checkpoint (`c0e526b`) on
  `codex/regression-provider-integration`.

### Task 2: Source-bound MEA records and immutable fit specification

**Files:**
- Modify: `src/epcsaft_regression/records.py`
- Modify: `src/epcsaft_regression/workflow.py`
- Modify: `src/epcsaft_regression/__init__.py`
- Test: `tests/test_records.py`

**Interfaces:**
- Consumes: the source identity, correlations, grid, 2B topology, bounds,
  scales, and starts frozen in `docs/science/general-parameter-regression.md`.
- Produces: `load_pure_saturation_dataset("monoethanolamine")` and the existing
  `PureSaturationFitSpecification` shape with five ordered parameter entries.

- [x] Add failing record tests for the exact PDF hash, DOI, 15 temperatures,
  calculated pressures, DIPPR-105 densities, molar-mass conversion, all-training
  partition, five parameter identities, and two declared starts.
- [x] Run the focused record tests and retain the expected unsupported
  component/specification failure.
- [x] Add the immutable MEA dataset/specification through the existing closed
  record owners; do not add a registry or external runtime data dependency.
- [x] Run all record tests and verify methane/ethane/propane serialized parity.
- [x] Commit the green source-contract checkpoint (`cd678b8`).

### Task 3: Data-sized existing Ceres owner

**Files:**
- Modify: `src/epcsaft_regression/native/pure_saturation_fit_internal.hpp`
- Modify: `src/epcsaft_regression/native/pure_saturation_contract.cpp`
- Modify: `src/epcsaft_regression/native/pure_saturation_fit.cpp`
- Modify: `src/epcsaft_regression/native/module.cpp`
- Modify: `src/epcsaft_regression/workflow.py`
- Test: `tests/test_native_fit.py`

**Interfaces:**
- Consumes: three- or five-parameter payload and the installed Provider SDK
  tail from Task 1.
- Produces: the existing `fit_pure_saturation` and
  `PureSaturationFitResult`, with variable-length `parameters`, exact
  smooth-observable Jacobians, rank diagnostics, and a separate bounded
  scientific-comparison status.

- [x] Add failing native tests proving a `30 x 5` MEA residual/Jacobian,
  start and terminal directional derivatives, the smooth-objective cost
  identity, full and parameter rank 5,
  and unchanged methane/ethane/propane results.
- [x] Run the focused native tests and retain the fixed-dimension/callback
  failure.
- [x] Replace fixed row/parameter arrays with vectors inside the existing
  owner, dispatch the three-parameter or seven-coordinate Provider callback,
  and assemble `du/dz=-H_u^{-1}H_z` for the MEA observable objective.
- [x] Make `PureSaturationFitResult.parameters` a variable-length tuple while
  preserving the same result type and all existing fields.
- [ ] Run the final native and Python fit matrix, including unsupported callback,
  fingerprint, topology, incomplete-column, active-bound, and multistart
  failures.
- [ ] Commit the green engine checkpoint.

### Task 4: Installed-artifact reconstruction evidence

**Files:**
- Modify: `docs/science/general-parameter-regression.md`
- Create: one canonical source/evidence data owner only if the existing record
  serializer cannot retain the 15 calculated rows and result diagnostics.
- Test: `tests/test_native_fit.py`

**Interfaces:**
- Consumes: exact installed Provider and Regression wheels.
- Produces: one deterministic Baygi reconstruction result with paper-style
  pressure/density AADs and descriptive parameter differences.

- [ ] Build and install one Provider wheel and one Regression wheel in an
  isolated environment; record exact wheel/header hashes.
- [ ] Run primary and confirmation starts and retain solver, numerical,
  physical/workflow, rank/conditioning, active-bound, AAD, and descriptive
  parameter-comparison fields.
- [ ] Verify exact callback and residual directional derivatives, deterministic
  output, wheel contents/import/linkage, and negative-space exclusions.
- [ ] Run the complete isolated installed-wheel suites and repository cleanup
  audits.
- [ ] Update only the canonical science/status owners with the measured
  candidate outcome and commit the bounded evidence checkpoint.
