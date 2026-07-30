# ePC-SAFT Regression

This repository owns one strict pure-saturation Ceres workflow. Its methane and
ethane forms are accepted reproducible workflows. Checkpoint A extends that
same workflow to one local propane candidate from the approved Glos 2004
direct-experimental packet. Every Helmholtz value and derivative comes from the
installed `epcsaft.native_sdk.v1` EOS SDK.

The public workflow is:

```python
from epcsaft import Mixture, Parameters
from epcsaft_regression import (
    ETHANE_SATURATION_FIT_V1,
    fit_pure_saturation,
    load_pure_saturation_dataset,
)

model = Mixture(
    Parameters.from_catalog(
        "gross-2001-methane-ethane",
        components=("ethane",),
        version=1,
    )
)
result = fit_pure_saturation(
    model=model,
    dataset=load_pure_saturation_dataset("ethane"),
    specification=ETHANE_SATURATION_FIT_V1,
)
```

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
no added runtime seam. The candidate excludes binary interactions,
association, electrolytes, reactions, generic target families, and parameter
persistence. See [the scientific contract](docs/science/pure-saturation-regression.md)
and [candidate capability record](evidence/candidate-capability.yaml).

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
