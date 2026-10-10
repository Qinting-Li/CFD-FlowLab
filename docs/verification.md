# Verification

## Formula checks

The standalone C++ test verifies:

- consistency with the physical Euler flux for equal states;
- Tadmor's entropy-conservation identity;
- non-positive two-state entropy production;
- positive entropy-Hessian quadratic forms;
- finite fluxes for near-equal and low-pressure states.

The JavaScript oracle sampled 10,000 admissible state pairs.

| Metric | Result |
|---|---:|
| Maximum Tadmor identity defect | `7.327e-15` |
| Maximum same-state entropy production | `0` |
| Minimum inverse-Hessian quadratic form | `1.182e-2` |

## Solver checks

The PHengLEI structured solver was built with MSVC and Microsoft MPI. Roe and entropy-stable flux smoke tests both completed two iterations on the 79,520-cell Mach 8.1 sphere case. The entropy-stable run then completed 25,000 iterations and reported `Job Done`.

## Limitation

The base flux uses raw cell states, while dissipation uses Minmod-reconstructed states. This hybrid state pairing does not provide a proof of global entropy stability for the complete semi-discretization. The repository claims formula consistency and successful numerical execution only.
