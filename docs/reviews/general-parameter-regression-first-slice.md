# General Parameter Regression First-Slice Evidence

Status: `FIT_READY_AUTHORITY_NEUTRAL`

Date: 2026-07-27

## Implemented contract

The first general family fits one shared, dimensionless `k_ij` for a
caller-supplied neutral, nonassociating binary model and source-bound
fixed-composition VLE rows. Each training row owns two lifted log-volume
coordinates and four scaled residuals:

```text
(P_liquid - P_observed) / pressure_scale
(P_vapor - P_observed) / pressure_scale
(mu_1,liquid/RT - mu_1,vapor/RT) / chemical_potential_scale_1
(mu_2,liquid/RT - mu_2,vapor/RT) / chemical_potential_scale_2
```

The exact `4R x (1 + 2R)` Jacobian consumes the installed Provider's
row-major value/gradient/Hessian contract in `(n1,n2,V,k_ij)`. No third
derivatives, density roots, copied EOS equations, numerical production
derivatives, or Equilibrium dependency are present.

## Exact subjects

- Provider commits:
  `651b778de5b25199a9efe96ac67b89d2fd113a0a`,
  `30f896df4f63a8b12dbc8047923306364c416966`,
  `5d9065110dae1f2548fbf831e92bd5c3362d58d1`
- Provider tree: `9badefc22bbb5d6736d77fc5bcd8524c2a73dc71`
- Provider wheel:
  `/home/tnnrpolley21/Workspaces/Engineering/ePC-SAFT-project/artifacts/provider-general-parameter-sdk-v1/5d90651/epcsaft-0.1.0.dev0-cp313-cp313-linux_x86_64.whl`
- Provider wheel SHA-256:
  `f5c06ee9e5dbe29bcca881378355672057283a69e59e3ec8561c0d1378c60d14`
- Installed Provider header SHA-256:
  `18d74d761290154dad3a336ced0d25c0b2badeab7e9462de3485a99e619fae1a`
- Regression commits:
  `626b62df1d7ba43275a3931748ad60458325609c`,
  `dd45126095689d01110e1152ea27a90e4a82694e`,
  `8816c3cd53e1f1b779e92ca43f669cc70d01fc78`
- Regression tree: `08af610318df82a7ee4f286f46d06b9bfce4105e`
- Regression wheel:
  `/home/tnnrpolley21/Workspaces/Engineering/ePC-SAFT-project/artifacts/regression-general-parameter-v1/8816c3c/epcsaft_regression-0.1.0.dev0-cp313-cp313-linux_x86_64.whl`
- Regression wheel SHA-256:
  `f2a8f757eb56b7b1c21676245d40ced5f58b95b3ff233b9c6a233b9061365adc`

## Verification

The exact installed wheels passed:

```text
25 focused installed-artifact tests passed in 0.90 s
87 full installed-artifact tests passed; 1 campaign deselected in 226.29 s
1 explicit all-17-row campaign passed in 1.07 s
```

The explicit campaign consumes
`evidence/may-2015-methane-ethane-vle.csv`, SHA-256
`5cd1e74925a3c6504f5106dcf911f2cae2d6e99a5133fccc20454d8991bdbc7f`.
It reports:

```text
rows                         17
residuals                    68
variables                    35
Ceres termination            CONVERGENCE
full lifted rank             35
projected parameter rank     1
fitted k_ij                  -0.00843032298906253
active parameter bound       none
declared confirmation starts 2
native fit elapsed           0.083270 s
```

Independent scientific and code reviews caught and closed false convergence
reporting, missing Provider-domain preflight, unbounded fixed-field reads,
post-Ceres diagnostic loss, incomplete result provenance, Provider capability
superset handling, per-iteration Regression allocations, legacy native-model
loading, and accepted pure-component numerical parity.

## Interpretation

This evidence establishes a reusable, source-bound, exact-derivative
fixed-composition neutral-binary `k_ij` runtime. It does not establish
predictive validity, global uniqueness, uncertainty, catalog authority, or
support for `l_ij`, association, `k_hb_ij`, polar, dielectric, reactive, or
other parameter families. The May rows are all training data. The historical
pressure-closure miss remains immutable negative evidence and is not relaxed
or reinterpreted here.
