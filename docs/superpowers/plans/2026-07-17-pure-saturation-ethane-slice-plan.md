# Pure Saturation Ethane Slice Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the methane-only regression surface with one strict pure-saturation workflow that reproduces methane and adds a source-backed ethane fit through the existing provider parameter-Hessian entry.

**Architecture:** Records and workflow become component-neutral only where a second real component proves the abstraction. One Ceres engine continues to own the lifted-volume residuals and exact Jacobian; immutable methane and ethane specifications own all component-specific facts.

**Tech Stack:** Python 3.13, C++17, Ceres 2.2, CPython C API, CMake, scikit-build-core, pytest, installed `epcsaft.native_sdk.v1`.

## Global Constraints

- Read `AGENTS.md`, `CONTEXT.md`, `ARCHITECTURE.yaml`, the accepted methane design/receipt, and the canonical Slice 2 design before editing.
- Correct stale pending-authority wording before describing a new candidate; do not alter the accepted receipt.
- Preserve the accepted methane numerical result through the generalized API.
- Use only the existing provider `evaluate_pure_phase_parameters` entry; no provider source, private imports, or copied EOS equations.
- Keep one native target, one workflow, one result representation, and one residual/Jacobian owner.
- Delete methane-only public aliases and files in the same candidate; add no compatibility wrapper or registry.
- Use only audited real NIST rows. Equal weights are numerical normalization, not uncertainty.
- Keep predictive admission, uncertainty, provider persistence, binary fitting, association, electrolytes, and generic regression closed.

---

### Task 1: Reconcile accepted authority and freeze the generic public contract

**Files:**
- Modify: `AGENTS.md`
- Modify: `CONTEXT.md`
- Modify: `ARCHITECTURE.yaml`
- Modify: `README.md`
- Modify: `tests/test_records.py`
- Modify: `tests/test_native_fit.py`

**Interfaces:**
- Consumes: accepted receipt `promotion-0020-regression-methane-saturation-v1`.
- Produces: exact public names `load_pure_saturation_dataset`, `fit_pure_saturation`, and `PureSaturationFitResult`.

- [ ] **Step 1: Correct stale authority prose**

Change only statements that call promotion 0020 pending. State that the clean regression repository is authoritative for the exact reproducible methane workflow, while predictive and fitted-parameter admission remain absent.

- [ ] **Step 2: Add failing public-surface tests**

Require these root exports and the absence of old names:

```python
assert hasattr(epcsaft_regression, "load_pure_saturation_dataset")
assert hasattr(epcsaft_regression, "fit_pure_saturation")
assert hasattr(epcsaft_regression, "PureSaturationFitResult")
for retired in (
    "load_methane_dataset",
    "fit_methane_saturation",
    "MethaneFitResult",
    "METHANE_FIT_SPECIFICATION_V1",
):
    assert not hasattr(epcsaft_regression, retired)
```

- [ ] **Step 3: Run tests and verify RED**

Run: `python -m pytest tests/test_records.py tests/test_native_fit.py -q`

Expected: failures because the generic names do not exist.

- [ ] **Step 4: Commit the authority correction and RED contract**

```bash
git add AGENTS.md CONTEXT.md ARCHITECTURE.yaml README.md tests/test_records.py tests/test_native_fit.py
git commit -m "test: freeze pure saturation regression surface"
```

### Task 2: Generalize records and retain the audited ethane table

**Files:**
- Modify: `src/epcsaft_regression/records.py`
- Create: `src/epcsaft_regression/data/ethane_saturation.csv`
- Modify: `src/epcsaft_regression/__init__.py`
- Modify: `tests/test_records.py`

**Interfaces:**
- Consumes: exact retained lab source `data/reference/pure_component/saturation_properties/ethane/saturation_properties.csv`.
- Produces: `PureSaturationDataset`, `PureSaturationFitSpecification`, `METHANE_SATURATION_FIT_V1`, `ETHANE_SATURATION_FIT_V1`, and the closed loader.

- [ ] **Step 1: Add exact ethane source tests**

Freeze all ten rows from 100 through 280 K, source URL identity, units, packaged SHA-256, and this partition:

```python
assert dataset.training_temperatures_k == (140.0, 180.0, 220.0, 260.0)
assert dataset.held_out_temperatures_k == (120.0, 160.0, 200.0, 240.0)
assert dataset.stress_temperatures_k == (100.0, 280.0)
```

Require `load_pure_saturation_dataset("methane")` and `("ethane")`; reject aliases, case variants, unknown strings, and non-strings.

- [ ] **Step 2: Run record tests and verify RED**

Run: `python -m pytest tests/test_records.py -q`

Expected: failure because generic records and ethane data are absent.

- [ ] **Step 3: Rename record types in place**

Rename methane-specific dataclasses to `PureSaturationDataset` and `PureSaturationFitSpecification`. Add immutable `component_id`, expected fingerprint, molar mass, and explicit training/held-out/stress temperature tuples. Keep the existing field-native units and exact source identity checks.

Implement the closed loader directly:

```python
def load_pure_saturation_dataset(component_id: str) -> PureSaturationDataset:
    if type(component_id) is not str:
        raise TypeError("component_id must be an exact string")
    if component_id == "methane":
        return _load_dataset("methane_saturation.csv", METHANE_SOURCE_V1)
    if component_id == "ethane":
        return _load_dataset("ethane_saturation.csv", ETHANE_SOURCE_V1)
    raise ValueError("component_id must be 'methane' or 'ethane'")
```

Do not add a mutable mapping, plugin registry, file-path argument, or fallback dataset.

- [ ] **Step 4: Copy and hash the audited ethane CSV**

Copy exactly the retained decimal strings and LF-normalized rows. Record both the lab-source hash and packaged hash in the immutable source identity. Do not retain unrelated pure-component tables.

- [ ] **Step 5: Run record tests**

Run: `python -m pytest tests/test_records.py -q`

Expected: all record tests pass.

- [ ] **Step 6: Commit records and data**

```bash
git add src/epcsaft_regression/records.py src/epcsaft_regression/data/ethane_saturation.csv src/epcsaft_regression/__init__.py tests/test_records.py
git commit -m "feat: add strict ethane saturation records"
```

### Task 3: Generalize the one native fit engine

**Files:**
- Rename: `src/epcsaft_regression/native/methane_contract.cpp` -> `src/epcsaft_regression/native/pure_saturation_contract.cpp`
- Rename: `src/epcsaft_regression/native/methane_fit.cpp` -> `src/epcsaft_regression/native/pure_saturation_fit.cpp`
- Rename: `src/epcsaft_regression/native/methane_fit.hpp` -> `src/epcsaft_regression/native/pure_saturation_fit.hpp`
- Rename: `src/epcsaft_regression/native/methane_fit_internal.hpp` -> `src/epcsaft_regression/native/pure_saturation_fit_internal.hpp`
- Modify: `src/epcsaft_regression/native/module.cpp`
- Modify: `CMakeLists.txt`
- Modify: `tests/test_native_fit.py`

**Interfaces:**
- Consumes: one immutable component-specific payload and provider parameter Hessian.
- Produces: one native evaluator/solver for either admitted pure component.

- [ ] **Step 1: Parameterize the native Jacobian test**

Run the centered directional residual check once for methane and once for ethane. Each case must construct its exact model, dataset, specification, source fingerprint, and training variables. The expected Jacobian remains provider Hessian plus exact affine/log/density chain rules.

- [ ] **Step 2: Run the directional tests and verify RED**

Run: `python -m pytest tests/test_native_fit.py -k jacobian -q`

Expected: the ethane case fails because the native identity is methane-specific.

- [ ] **Step 3: Rename files and remove methane identity from the engine**

Use the generic file names above. Replace hard-coded methane identity fields with payload fields that are already validated by exact Python records: component ID, expected provider fingerprint, molar mass, row IDs, and row partition. Keep parameter coordinate names `(m,sigma,epsilon/k)` unchanged.

The C++ parser must reject a payload whose component ID, source identity, row count/order, fingerprint, units, or specification identity differs from the corresponding immutable Python specification. Do not infer a component from molar mass or fingerprint.

- [ ] **Step 4: Update the one native target**

Change only the three source paths in `CMakeLists.txt`. Do not add a second library or target.

- [ ] **Step 5: Run native tests**

Run: `python -m pytest tests/test_native_fit.py tests/test_transport.py -q`

Expected: methane and ethane derivative/transport tests pass; malformed identity tests fail loudly as expected.

- [ ] **Step 6: Commit the native generalization**

```bash
git add CMakeLists.txt src/epcsaft_regression/native tests/test_native_fit.py tests/test_transport.py
git commit -m "refactor: generalize pure saturation fit engine"
```

### Task 4: Generalize the Python workflow and prove methane parity

**Files:**
- Modify: `src/epcsaft_regression/workflow.py`
- Modify: `src/epcsaft_regression/__init__.py`
- Modify: `tests/test_native_fit.py`
- Modify: `tools/run_candidate.py`

**Interfaces:**
- Consumes: generic records and native solver.
- Produces: `fit_pure_saturation(...) -> PureSaturationFitResult` and reproducible component-specific candidate receipts.

- [ ] **Step 1: Freeze methane parity before renaming**

Record the accepted result values from `evidence/candidate-fit-receipt.json`: final parameters, initial/final cost, ranks, and all nine reporting predictions. Compare with existing accepted tolerances; do not regenerate expected values after implementation.

- [ ] **Step 2: Rename workflow types and function**

Rename `MethaneFitResult` to `PureSaturationFitResult` and `fit_methane_saturation` to `fit_pure_saturation`. Replace exact methane type checks with exact generic record checks and verify `dataset.component_id`, `specification.component_id`, model fingerprint, and source identity all match before calling native code.

- [ ] **Step 3: Generalize the candidate runner**

Add required `--component methane|ethane`; bind the selected source/specification, artifact hashes, and output receipt subject. The runner writes one component-specific receipt and never loops over an implicit catalog.

- [ ] **Step 4: Run methane parity and ethane smoke tests**

Run: `python -m pytest tests/test_native_fit.py -q`

Expected: accepted methane values remain within frozen tolerances; ethane returns finite bounded parameters, full parameter rank, usable confirmation, and valid reporting states.

- [ ] **Step 5: Commit the public workflow**

```bash
git add src/epcsaft_regression/workflow.py src/epcsaft_regression/__init__.py tests/test_native_fit.py tools/run_candidate.py
git commit -m "feat: expose pure saturation regression workflow"
```

### Task 5: Installed-artifact real-data candidate and documentation

**Files:**
- Modify: `docs/science/methane-saturation-regression.md`
- Create: `docs/science/pure-saturation-regression.md`
- Modify: `README.md`
- Modify: `CONTEXT.md`
- Modify: `ARCHITECTURE.yaml`
- Modify: `evidence/candidate-capability.yaml`
- Create: `evidence/ethane-candidate-fit-receipt.json`

**Interfaces:**
- Consumes: immutable provider and regression wheels.
- Produces: one ethane candidate evidence subject without predictive or persistence admission.

- [ ] **Step 1: Replace methane-only science ownership**

Move shared equations, units, transforms, residuals, and derivative contract into `pure-saturation-regression.md`. Keep methane and ethane source partitions explicit. Remove the superseded methane-only science file after all live links are updated.

- [ ] **Step 2: Build the reproducible wheel**

Build with fixed `SOURCE_DATE_EPOCH`, install the exact provider and regression wheels in an isolated Python 3.13 environment, verify every installed non-RECORD member against its wheel, and run:

```text
python -m pytest -q
python tools/run_candidate.py --component methane ...
python tools/run_candidate.py --component ethane ...
```

Expected: all tests pass, methane parity passes, and ethane produces a deterministic candidate receipt.

- [ ] **Step 3: Record honest evidence boundaries**

Report training, held-out, and stress metrics separately. Do not set physical accuracy success from Ceres success, objective decrease, or stress endpoints. Keep provider persistence, uncertainty, global identifiability, and generic-family claims false.

- [ ] **Step 4: Run final ratchet and cleanup checks**

Require one native target and no retired methane symbols/files, aliases, registries, private provider imports, copied equations, sibling paths, or extra generated output. Run `git diff --check`, full tests, exact tracked-file review, and the cleanup hook.

- [ ] **Step 5: Commit and stop for independent review**

```bash
git add README.md CONTEXT.md ARCHITECTURE.yaml docs/science evidence
git commit -m "docs: record ethane saturation regression candidate"
```

Report commit/tree/wheel/receipt hashes and commands. Do not push, persist fitted parameters, broaden capability, or accept authority.
