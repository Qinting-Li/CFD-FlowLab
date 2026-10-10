# Entropy-Stable Flux for Compressible Euler Flow

This project studies an entropy-conservative Euler flux with entropy-variable dissipation and its integration into a structured finite-volume solver.

## Method

```text
f_num = f_EC(qL0, qR0)
        - 0.5 * lambda_max * Hbar_inverse * (vR_rec - vL_rec)
```

- `f_EC` uses logarithmic means of the raw cell states.
- Dissipation uses Minmod-reconstructed entropy variables.
- `lambda_max` is the maximum local acoustic wave speed.
- `Hbar_inverse` is the analytic inverse entropy Hessian.
- Roe flux is retained as the reference scheme.

## Verification

| Check | Result |
|---|---:|
| Maximum Tadmor identity defect | `7.327e-15` |
| Maximum sampled same-state entropy production | `0` |
| Minimum inverse-Hessian quadratic form | `1.182e-2` |
| C++ formula test | Passed |
| Mach 8.1 sphere run | 25,000 iterations, `Job Done` |

The hybrid use of raw and reconstructed states does not prove global entropy stability for the complete semi-discretization. See [docs/verification.md](docs/verification.md).

## Run the standalone checks

```bash
cmake -S . -B build
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
node tests/check_entropy_stable_euler.js
```

## Repository layout

- `include`: standalone flux implementation
- `tests`: C++ and JavaScript verification
- `patches`: PHengLEI integration patch
- `examples/sphere`: minimal case, grid, results, and wall-pressure visualization

## Integration

Apply `patches/phenglei-integration.patch` to the compatible PHengLEI source tree, build PHengLEI, and run the solver from `examples/sphere`.

PHengLEI is an external host framework. This repository does not relicense or redistribute its full source tree.
