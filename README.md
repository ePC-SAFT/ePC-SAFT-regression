# ePC-SAFT Regression

This repository owns source-bound ePC-SAFT parameter fitting with Ceres. Every
thermodynamic primitive and EOS-explicit production derivative comes from an
immutable installed EOS artifact through `epcsaft.native_sdk.v1`; downstream
composed derivatives may additionally consume exact implicit sensitivities from
an installed Equilibrium artifact or another authorized evaluator. Regression
owns target contracts, parameter transforms, residual/Jacobian assembly, fit
diagnostics, and authority-neutral results without reimplementing those owners.

The current public surface includes ten scalar parameter families, a joint pure
`(m, sigma, epsilon/k)` block, a descriptor-driven fixed-topology association
block, and a specialized transport for downstream-composed positive
observations. Association coordinates are typed by component/site endpoints;
the installed EOS owns the finite topology, slot order, multiplicities,
physical domain, artifact identity, and exact derivative contract.

Production authority is narrower than executable mechanics. Only the exact
methane and ethane pure-saturation workflows named below are accepted
reproducible workflows. The general families and fixed-topology association
block are authority-neutral; their mechanics evidence does not identify or
admit a unique association topology or parameter tuple.

The accepted pure-saturation presentation workflow remains:

```python
from epcsaft import Mixture, Parameters
from epcsaft_regression import (
    ETHANE_SATURATION_FIT_V1,
    fit_pure_saturation,
    load_pure_saturation_dataset,
)

model = Mixture(
    Parameters.from_bundle(
        "/path/to/ePC-SAFT-data/packets/gross-2001-methane-ethane/1/parameters",
        components=("ethane",),
    )
)
result = fit_pure_saturation(
    model=model,
    dataset=load_pure_saturation_dataset("ethane"),
    specification=ETHANE_SATURATION_FIT_V1,
)
```

Caller-supplied datasets use the public
[`prepare -> preflight -> fit -> export` quickstart](docs/general-regression-quickstart.md).
It covers strict records, acquisition/objective provenance, grouped controls,
installed-capability resolution, local-rank preflight, deterministic export,
and every direct-EOS observation template. Direct `RegressionProblem` and
`fit_parameters` construction remains available with the same native payload.

Parameter inputs are explicit. Source-traceable reusable bundles live in
ePC-SAFT-data and are loaded with `Parameters.from_bundle(path)`; compact
caller-owned numerical inputs use `Parameters.from_dictionary(mapping)`.
Regression never locates a sibling Data checkout or reconstructs an installed
catalog.

Builds require Ceres 2.2 and `EPCSAFT_INCLUDE_DIR` set to the public include
directory of an installed compatible EOS wheel. The build does not search
sibling source trees.

Accepted historical receipts `promotion-0020-regression-methane-saturation-v1`
and `promotion-0023-regression-pure-saturation-ethane-v1` make this repository
the production owner of the exact reproducible methane and ethane workflows;
`state-0025-regression-ethane-publication` verifies ethane publication. They do
not admit fitted parameters as predictive or scientific authority. Held-out
errors are descriptive because no admission cutoff was approved, and the
reporting-block directional Jacobian remains an explicit evidence limit with
no added runtime seam. Those receipts do not promote the authority-neutral
general parameter families, fixed-topology association mechanics, downstream-composed evaluator
transport, or any fitted parameter bundle. Regression still excludes reactive
chemistry/equilibrium algorithms, polar regression, arbitrary residual plugins,
EOS catalog persistence, global-identifiability claims, and predictive
authority without separate evidence.

See the
[general parameter-family contract](docs/science/general-parameter-regression.md),
[literature reproduction contract](docs/science/literature-reproduction-contract.md),
[accepted pure-saturation contract](docs/science/pure-saturation-regression.md),
[MEA coordination contract](docs/science/mea-coupled-regression-master-plan.md),
and [candidate capability record](evidence/candidate-capability.yaml).

Regression's routine suite keeps compact numerical-tolerance sentinels for the
main scalar, pair, joint-pure, fixed-topology association, and composed-observation workflows.
Run it with `.venv/bin/pytest -q`; the checked-in pytest configuration excludes
the explicit `campaign` marker.
Full literature reproductions, broad benchmark matrices, identifiability and
uncertainty campaigns, and induced/cross-association comparisons belong as
durable installed-artifact evidence in the sibling Validation repository.

Project doctrine and cross-repository coordination are owned by
[ePC-SAFT Management](https://github.com/ePC-SAFT/ePC-SAFT-management) under
doctrine revision 5. EOS is the canonical active role name; retained
`provider-*` artifact paths and receipt identifiers remain immutable
historical evidence.

Codex-created worktrees run the checked-in setup at
`.codex/environments/setup.sh`. It verifies and installs the exact retained EOS
wheel before building Regression; `EPCSAFT_EOS_WHEEL` may point to that same
hash-bound artifact when the project artifact root is elsewhere.

The local propane Checkpoint A candidate converges and confirms with complete
rank, but held-out 120 K pressure closure is about `3.29e-8`, above the frozen
`1e-8` numerical gate. It is therefore retained as blocked evidence, not as a
predictive or promoted workflow. The pointwise Glos uncertainties remain
reporting context and are not substituted for that gate.
