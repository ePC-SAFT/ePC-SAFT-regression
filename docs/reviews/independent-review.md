# Independent Review: First Methane Regression Slice

Reviewer: `/root/independent_review`
Review date: 2026-07-17
Authority: implementation-quality review only; no scientific admission or runtime promotion

## Initial verdict

The first review blocked commit because the numerical result was promising but
the provenance and acceptance contracts could describe a different problem
than the one executed. It verified the provider/regression ownership split,
the Hessian/Jacobian equations, source hashes, isolated tests, wheel contents,
and absence of copied provider symbols.

The blocking findings were:

1. supplied wheel hashes were not bound to installed runtime files;
2. the compiled native problem dropped source and row identity;
3. reporting transforms, bounds, confirmation starts, and solver controls were hidden constants;
4. source records accepted changed metadata and row-ID/temperature swaps;
5. reporting physical validity lacked an equilibrium-closure threshold;
6. the required architecture baseline was absent;
7. capsule-prefix inspection used an incompatible C++ object type.

Additional findings requested structural Jacobian completeness, portable and
more comprehensive receipt hashing, explicit phase densities, source notes for
the retained start/molar mass, and splitting the 875-line native owner.

## Adjudication

All findings were accepted and addressed:

- The candidate runner hashes artifacts before package import, byte-compares
  every non-`RECORD` wheel member to the installed distribution, then checks
  import origins. The receipt records portable member names and a canonical
  hash over its full payload.
- The native contract independently retains and exactly validates source
  metadata/hashes, row IDs and values, units, names, transforms, bounds, fixed
  amount, provider commit/wheel/fingerprint, and solver controls. Native row IDs
  and the full compiled identity round-trip to Python and are checked there.
- Every reporting, confirmation, acceptance, and Ceres control is an immutable
  specification field included in the native payload and receipt.
- Exact source fields and row-ID/temperature mappings are rejected on mutation.
- Reporting pressure and chemical-potential closures have explicit `1e-8`
  numerical gates, separate from predictive accuracy.
- Capsule prefix fields are copied with `memcpy` before the tail size gate.
- Every Jacobian entry is explicitly initialized through the assembly path;
  rank and condition diagnostics remain separate from structural completeness.
- Training and reporting records include both phase molar and mass densities.
- The native contract parser was separated from the residual/Ceres owner while
  retaining one extension target and no additional public surface.

## Follow-up verification

The reviewer inspected the live fixes, independently reproduced 17 installed-
wheel tests, the receipt, its canonical digest, artifact binding, architecture
counts, and negative-space scans, and cleared every substantive finding. The
user-authorized `.idea/` project configuration is retained unchanged; its
existing local `workspace.xml` ignore remains in force.

Final readiness: **READY** for local commit and candidate manager review only; no scientific admission or production authority is granted.
