# ePC-SAFT Regression

This repository owns one strict pure-saturation Ceres workflow. Its methane and
ethane forms are accepted reproducible workflows. Checkpoint A extends that
same workflow to one local propane candidate from the approved Glos 2004
direct-experimental packet. Every Helmholtz value and derivative comes from the
installed `epcsaft.native_sdk.v1` provider capsule.

The public workflow is:

```python
from epcsaft import EPCSAFT, ParameterBundle
from epcsaft_regression import (
    ETHANE_SATURATION_FIT_V1,
    fit_pure_saturation,
    load_pure_saturation_dataset,
)

model = EPCSAFT(
    ParameterBundle.from_catalog("gross-2001-methane-ethane", version=1).select(
        ("ethane",)
    )
)
result = fit_pure_saturation(
    model=model,
    dataset=load_pure_saturation_dataset("ethane"),
    specification=ETHANE_SATURATION_FIT_V1,
)
```

Builds require Ceres 2.2 and `EPCSAFT_INCLUDE_DIR` set to the public include
directory of an installed compatible provider wheel. The build does not search
sibling source trees.

Accepted migration receipts `promotion-0020-regression-methane-saturation-v1`
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
