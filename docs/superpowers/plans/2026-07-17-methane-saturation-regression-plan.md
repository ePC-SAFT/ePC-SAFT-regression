# Pure-Methane Saturation Regression Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build and verify the first clean, installed-provider-only Ceres regression candidate for methane `m`, `sigma`, and `epsilon/k`.

**Architecture:** Python owns strict immutable source/problem records and the single public workflow. One native CPython module owns lifted-volume Ceres residuals, exact provider-Hessian Jacobians, solve acceptance, diagnostics, and fixed-parameter reporting solves through the model-bound capsule.

**Tech Stack:** Python 3.13, C++17, CPython C API, CMake, scikit-build-core, Ceres 2.2, Eigen through Ceres, installed `epcsaft` wheel and public `native_sdk_v1.h`.

## Global Constraints

- Use provider commit `4b10cb899c94687cae734980285badb224dc95e6` and wheel SHA-256 `f92f79c8d6f614660e5c201b7061c9b02b5cd1a25a4ed8c8fee0b59adaabf2bf`.
- Fit exactly methane `m`, `sigma [angstrom]`, and `epsilon/k [K]` at 110, 130, 150, and 170 K with fixed `n = 1 mol`.
- Use only `evaluate_pure_phase_parameters` and exact Hessian/linear chain-rule Jacobians in production.
- No provider source path, private provider header/import, copied equation, duplicate provider symbol, alternate derivative backend, persistence, or broader regression family.
- Separate solver convergence, numerical convergence, and physical validity; candidate evidence does not self-admit capability.

---

### Task 1: Strict source and problem records

**Files:**
- Create: `src/epcsaft_regression/records.py`
- Create: `src/epcsaft_regression/data/methane_saturation.csv`
- Create: `tests/test_records.py`

**Interfaces:**
- Produces: immutable `SourceIdentity`, `SaturationObservation`, `MethaneSaturationDataset`, `MethaneFitSpecification`, and `load_methane_dataset()`.

- [ ] **Step 1: Write failing contract tests.** Assert the exact nine source rows, exact training IDs, units/source/hash, tuple ordering, explicit starts/bounds/scales, and rejection of missing units, duplicate IDs, wrong order, nonfinite values, changed source identity, and unsupported species.
- [ ] **Step 2: Run RED.** Run `python3.13 -m pytest tests/test_records.py -q`; expect import failure because `epcsaft_regression.records` does not exist.
- [ ] **Step 3: Implement the immutable records and exact CSV.** Parse the fixed header with `csv.reader`, construct frozen slotted dataclasses, validate every invariant in `__post_init__`, and return tuple-backed records only.
- [ ] **Step 4: Run GREEN.** Run `python3.13 -m pytest tests/test_records.py -q`; expect all record tests to pass.

### Task 2: Package and capsule boundary

**Files:**
- Create: `pyproject.toml`
- Create: `CMakeLists.txt`
- Create: `src/epcsaft_regression/native/methane_fit.hpp`
- Create: `src/epcsaft_regression/native/module.cpp`
- Create: `tests/test_transport.py`

**Interfaces:**
- Consumes: the installed public header `epcsaft/native_sdk_v1.h` from explicit `EPCSAFT_INCLUDE_DIR`.
- Produces: private module `epcsaft_regression._native` with `transport_info(capsule)` and later `evaluate`/`solve` functions.

- [ ] **Step 1: Write the failing cross-wheel test.** In an environment containing the pinned provider wheel, construct the exact methane `EPCSAFT` model and assert the extension validates ABI version, minimum table size, parameterized result size, and a non-null parameterized evaluator.
- [ ] **Step 2: Run RED.** Configure/build the extension and run the focused test; expect failure because the native target/module is absent.
- [ ] **Step 3: Implement the minimum package shape and ABI gate.** Require the explicit installed-header directory, compile one module target, use `PyCapsule_GetPointer`, copy only the prefix before the size gate, then validate the parameterized tail and retain the capsule during calls.
- [ ] **Step 4: Run GREEN.** Rebuild and run `python -m pytest tests/test_transport.py -q`; expect the ABI gate to pass against the pinned installed wheel.

### Task 3: Exact lifted residual/Jacobian evaluator

**Files:**
- Create: `src/epcsaft_regression/native/methane_fit.cpp`
- Modify: `src/epcsaft_regression/native/module.cpp`
- Create: `tests/test_native_fit.py`

**Interfaces:**
- Consumes: ordered training rows and the immutable fit specification.
- Produces: private native `evaluate(capsule, payload, variables)` returning residuals, row diagnostics, and the complete row-major Jacobian.

- [ ] **Step 1: Write failing anchor and directional tests.** Freeze provider-direct values at the 130 K start state, assert raw/scaled residual ordering, then compare centered residual directional differences against `Jv` for a fixed non-axis direction.
- [ ] **Step 2: Run RED.** Run the focused native test; expect missing native `evaluate`.
- [ ] **Step 3: Implement the evaluator.** Map affine parameter and log-volume coordinates, call the provider twice per row, assemble the four scaled residuals and all exact Jacobian columns, and reject provider errors, nonfinite values, nonpositive `Phi_VV`, or lost phase ordering.
- [ ] **Step 4: Run GREEN.** Rebuild and run the focused anchor/directional tests; require the independent anchor tolerances and directional error thresholds to pass.

### Task 4: Ceres fit, diagnostics, and reporting solve

**Files:**
- Modify: `src/epcsaft_regression/native/methane_fit.cpp`
- Modify: `src/epcsaft_regression/native/module.cpp`
- Create: `src/epcsaft_regression/workflow.py`
- Create: `src/epcsaft_regression/__init__.py`
- Modify: `tests/test_native_fit.py`

**Interfaces:**
- Produces: `fit_methane_saturation(model, dataset, specification) -> MethaneFitResult` with immutable diagnostics and reporting rows.

- [ ] **Step 1: Write the failing workflow test.** Require deliberate parameter movement, respected bounds, complete Jacobian columns, singular/rank diagnostics, per-row residuals, separate convergence statuses, stable distinct phases, all nine reporting rows, and explicit failure reasons.
- [ ] **Step 2: Run RED.** Run the focused workflow test; expect missing public function/result.
- [ ] **Step 3: Implement the minimal Ceres workflow.** Add the 16-residual/11-variable exact cost, physical bounds, convergence gate, perturbed-volume confirmation solve, Eigen singular diagnostics, and nine fixed-parameter three-variable saturation reporting solves. Convert the private tuple transport into one frozen public result representation.
- [ ] **Step 4: Run GREEN.** Rebuild and run all three tests; require the candidate fit and every structural/physical assertion to pass.

### Task 5: Scientific record, candidate receipt, and architecture baseline

**Files:**
- Create: `docs/science/methane-saturation-regression.md`
- Create: `tools/run_candidate.py`
- Create: `evidence/candidate-capability.yaml`
- Create: `evidence/candidate-fit-receipt.json`
- Create: `ARCHITECTURE.yaml`
- Modify: `README.md`
- Modify: `CONTEXT.md`

**Interfaces:**
- Consumes: the installed provider and completed public workflow.
- Produces: a replayable candidate receipt and non-authoritative capability/architecture records.

- [ ] **Step 1: Write the candidate runner.** Require explicit `--provider-wheel`, verify its SHA-256 before import, execute the fit, and serialize exact provider/source/artifact identities, tuples, transforms, bounds, rows, diagnostics, held-out metrics, exclusions, and `candidate` authority status.
- [ ] **Step 2: Execute the candidate.** Run from an isolated environment containing only the pinned provider and regression artifacts; retain the JSON output and exact command/environment fields.
- [ ] **Step 3: Write the scientific and architecture records.** Document the equations consumed from the provider contract, units, transforms, residual scales/weights, derivative identities, validation boundary, public exports, dependencies, native target, line/file counts, and excluded scope.
- [ ] **Step 4: Validate records.** Parse YAML/JSON, verify referenced hashes and row coverage, and confirm no record claims accepted production or validation authority.

### Task 6: Source and isolated-wheel verification

**Files:**
- Modify only files implicated by verification failures.

**Interfaces:**
- Produces: final local evidence and a clean commit on `main`.

- [ ] **Step 1: Run source verification.** Build against a provider-only virtual environment and run the complete test suite, `compileall`, CMake build, compiler warnings, `git diff --check`, and static forbidden-pattern/symbol scans.
- [ ] **Step 2: Run isolated-wheel verification.** Build the regression wheel, create a fresh environment, install the pinned provider wheel and regression wheel only, rerun the tests and candidate runner, and hash both wheels plus the receipt.
- [ ] **Step 3: Run independent review.** Give a reviewer the exact diff, design, provider contract, test output, receipt, and negative-space scans; adjudicate every actionable finding and rerun affected checks.
- [ ] **Step 4: Run cleanup and commit.** Execute `bash "$HOME/.codex/hooks/codex-cleanup.sh" --repo-root .`, preserve pre-existing `.serena/`, verify the intended clean status, and commit locally on `main` without pushing.

## Self-review

- Spec coverage: every user requirement maps to Tasks 1–6; black-box scientific admission and catalog persistence remain explicitly excluded.
- Placeholder scan: no TBD/TODO or undefined follow-on implementation remains.
- Type consistency: the public workflow consumes the Task 1 dataset/specification records and returns the single Task 4 immutable result; private native tuples are transport only.
- Scope: one pure methane family, one native target, one provider dependency, and one candidate receipt.
