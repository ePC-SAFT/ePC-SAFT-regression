# ePC-SAFT Regression

This repository owns one strict pure-saturation Ceres workflow. It reproduces
the accepted methane fit and adds a local ethane candidate from retained NIST
saturation rows. It calls the installed `epcsaft.native_sdk.v1` provider
capsule for every Helmholtz value and derivative.

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

Accepted migration receipt `promotion-0020-regression-methane-saturation-v1`
makes this repository the production owner of the exact reproducible methane
workflow. It does not admit the fitted parameters as predictive or scientific
authority. Held-out errors are descriptive because no admission cutoff was
approved, and the
reporting-block directional Jacobian remains an explicit evidence limit with
no added runtime seam. The candidate excludes binary interactions,
association, electrolytes, reactions, generic target families, and parameter
persistence. See [the scientific contract](docs/science/pure-saturation-regression.md)
and [candidate capability record](evidence/candidate-capability.yaml).
