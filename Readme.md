# Entropy-Stable Numerical Flux for Compressible CFD

This repository studies the design, implementation, and verification of an entropy-stable numerical flux for compressible Euler flows. PHengLEI is used only as the host CFD framework for integration and flow-case evaluation.

## Research Core

The numerical flux is

```text
f_num = f_EC(qL^0, qR^0)
        - 0.5 * lambda_max * Hbar^{-1} * [v(qR^rec) - v(qL^rec)].
```

- `f_EC` is an entropy-conservative two-point flux based on logarithmic means.
- `qL^0` and `qR^0` are raw cell-centred states used by the base flux.
- `qL^rec` and `qR^rec` are Minmod-reconstructed states used by the entropy-variable jump and wave speed.
- `Hbar^{-1}` is the analytic inverse entropy Hessian at the averaged reconstructed state.
- Roe flux remains available as the reference scheme.

## Verification

| Check | Result |
|---|---:|
| Maximum Tadmor identity defect over 10,000 state pairs | `7.33e-15` |
| Maximum sampled same-state entropy production | `0` |
| Minimum sampled entropy-Hessian quadratic form | `1.18e-2` |
| C++ formula test | `1/1 passed` |
| Mach 8.1 sphere smoke tests | Roe and proposed flux completed |
| Full sphere run | `25,000` iterations, `Job Done` |

The full case uses a 79,520-cell structured sphere grid. Tecplot data, residual histories, surface quantities, and an SVG wall-pressure visualization are provided in `PHengLEI/examples/results`.

## Scope and Limitation

The entropy-conservative identity and local dissipation checks pass to numerical precision. Because the base flux uses raw states while dissipation uses reconstructed states, global entropy stability of the complete reconstructed semi-discretisation is not claimed. The numerical evidence and this limitation are documented in `API/Physics/tests/EntropyStableEulerReport.md`.

## Key Files

- `API/Physics/src/Flux_Inviscid.cpp`: flux implementation
- `CFD/NSSolver/src/NSSolverStruct.cpp`: structured-solver integration
- `API/Physics/tests`: analytical and numerical verification
- `PHengLEI/examples`: sphere case, outputs, and visualization

Run the independent numerical checks with:

```bash
node API/Physics/tests/EntropyStableEulerNumericalChecks.js
```

## Attribution

The method implementation and experiments are built on the PHengLEI CFD framework. The original framework and its licence remain acknowledged in this repository.
