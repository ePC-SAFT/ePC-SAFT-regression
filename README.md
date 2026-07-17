# ePC-SAFT Regression

This repository contains one candidate regression slice: a Ceres fit of pure
methane segment count, segment diameter, and dispersion energy to four retained
NIST saturation rows. It calls the installed `epcsaft.native_sdk.v1` provider
capsule for every Helmholtz value and derivative.

The public workflow is:

```python
from epcsaft import EPCSAFT, ParameterBundle
from epcsaft_regression import (
    METHANE_FIT_SPECIFICATION_V1,
    fit_methane_saturation,
    load_methane_dataset,
)

model = EPCSAFT(
    ParameterBundle.from_catalog("gross-2001-methane-ethane", version=1).select(
        ("methane",)
    )
)
result = fit_methane_saturation(
    model=model,
    dataset=load_methane_dataset(),
    specification=METHANE_FIT_SPECIFICATION_V1,
)
```

Builds require Ceres 2.2 and `EPCSAFT_INCLUDE_DIR` set to the public include
directory of an installed compatible provider wheel. The build does not search
sibling source trees.

This candidate has not received installed-artifact predictive admission and
does not change runtime authority. It excludes binary interactions,
association, electrolytes, reactions, generic target families, and parameter
persistence. See [the scientific contract](docs/science/methane-saturation-regression.md)
and [candidate capability record](evidence/candidate-capability.yaml).
