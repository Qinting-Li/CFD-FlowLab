# Entropy-stable Euler flux: implementation and verification

## Scope and existing implementation

The structured solver already reconstructs the primitive face states
`qL`/`qR` in `NSSolverStruct::GetInvFaceVar`, using the existing Minmod
limiter. The Roe flux remains available and remains the default. The new
`entropy_stable` selector is opt-in and requires that existing structured
Minmod reconstruction.

The original Roe entropy-fix adjusts Roe eigenvalues; it is not an
entropy-conservative flux. This change adds a separate ideal-gas,
single-component Euler entropy-conservative base flux and the requested
entropy-scaled dissipation. It does not change the Minmod algorithm or time
integration.

## Formula and state separation

For each face, the implementation keeps two distinct pairs of states:

* `qL^0`, `qR^0`: raw cell-centered primitive states, used only by the
  entropy-conservative base flux.
* `qL^rec`, `qR^rec`: existing Minmod-reconstructed primitive states, used
  for the entropy-variable jump and the maximum wave speed.

The implemented flux is

```text
f_num = f_EC(qL^0, qR^0; n)
        - 0.5 lambda_max Hbar^{-1}
          [v(qR^rec) - v(qL^rec)].
```

Here `Hbar^{-1}` is the analytic inverse entropy Hessian evaluated at the
arithmetic average of the reconstructed primitive states. It is the
entropy-scaled matrix product in the requested expression
(`Rtilde Rtilde^T = Hbar^{-1}`). `lambda_max` is the maximum of
`|u.n| + c` on the two reconstructed states. The flux normalizes `n`; the
solver applies face area later.

The entropy variables are formed from each reconstructed physical state,
not from a conservative-variable difference:

```text
s = ln(p) - gamma ln(rho)
beta = rho / (2 p)
v = [(gamma-s)/(gamma-1) - beta |u|^2,
     2 beta u, 2 beta v, 2 beta w, -2 beta].
```

The entropy-conservative Euler flux uses logarithmic means and satisfies
Tadmor's two-state identity to roundoff in the numerical checks.

The implementation validates the supported physical configuration (static
mesh, no preconditioning, five Euler equations, one species and one
temperature, matching ideal-gas `gamma`) and positive density and pressure.
The feature is intended for the structured solver; it is not a replacement
for the unstructured solver's separate limiter path.

## Applying the PHengLEI integration

This GitHub branch is based on a repository whose `main` contains no PHengLEI
source tree. The standalone header and tests can be reviewed here; the solver
integration is provided as
[`entropy-stable-euler-integration.patch`](./entropy-stable-euler-integration.patch).
Apply that patch from the root of a complete PHengLEI checkout:

```powershell
git apply API\Physics\tests\entropy-stable-euler-integration.patch
```

The patch changes only the scheme identifier, structured inviscid-flux
declaration and implementation, and structured solver dispatch/state access.
It was generated against the PHengLEI source baseline used for the local
build; this small GitHub repository cannot build the complete solver by
itself.

## Entropy-stability qualification

For a flux whose base and dissipative terms use the **same** two states, the
entropy-conservative identity and positive definiteness of `Hbar^{-1}` give
the usual two-state entropy inequality. The requested hybrid here deliberately
uses raw states for `f_EC` and reconstructed states for dissipation. Its
cell-center entropy production contains the cross term

```text
-0.5 lambda_max
 [v(qR^0)-v(qL^0)]^T Hbar^{-1}
 [v(qR^rec)-v(qL^rec)],
```

which is not the nonnegative quadratic form obtained when both jumps are
identical. Componentwise Minmod on primitive variables does not guarantee
that the reconstructed entropy jump has the sign of the raw entropy jump.
The deterministic scan below confirms that componentwise sign mismatches
occur. A separate scan found no positive cross entropy production in its
sampled states, but sampling is not a proof. Therefore this implementation
is **not claimed to be entropy stable as a complete reconstructed
semidiscretization**. No clipping or alternate reconstruction was introduced;
that would change the requested flux or the existing limiter.

## Build and automated checks

The complete Release project built successfully with MSVC 19.50 / Visual
Studio 2026 and Microsoft MPI. The executable was produced at:

```text
%TEMP%\phenglei-full-build\PHengLEIv3d0\Release\PHengLEIv3d0.exe
```

The independent C++ formula test executable completed successfully (CTest:
1/1 passed; direct run also passed). Its reported values were:

| C++ formula-test output | Result |
|---|---:|
| Tadmor identity defect | `-1.11022e-16` |
| Two-state entropy production | `-1.06766` |
| Reconstructed entropy dissipation quadratic form | `1.10064` |
| Minimum tested inverse-Hessian quadratic form | `0.0146709` |

The deterministic JavaScript numerical oracle can be rerun from the
repository root:

```powershell
node API\Physics\tests\EntropyStableEulerNumericalChecks.js
```

Observed results:

| Check | Result |
|---|---:|
| Random admissible two-state pairs | 10,000 |
| Maximum Tadmor identity defect | `7.327471962526033e-15` |
| Maximum entropy production when base and dissipation use the same pair | `0` |
| Minimum sampled `delta-v^T Hbar^-1 delta-v` | `0.01182242505351327` |
| Primitive-Minmod entropy-jump component sign mismatches | `23,310 / 500,000` |
| Hybrid raw-state entropy-production violations in random scan | `0 / 100,000` |
| Maximum sampled hybrid raw-state entropy production | `-0.006203814903235572` |

The last two rows are finite-sample observations only and do not establish a
general entropy inequality.

## Actual PHengLEI runs

The repository's current 3dball example was copied to two isolated temporary
directories. Each used the same structured grid, Mach 8.1, laminar physical
model, existing Minmod reconstruction, and all other copied example
parameters. Only the inviscid selector differed (`roe` versus
`entropy_stable`). In the temporary configuration copies only,
`maxSimuStep` was reduced to 2 and `intervalStepRes` to 1. The source example
files were not changed for these runs.

Reproduction outline from the repository root:

```powershell
$exe = Join-Path $env:TEMP 'phenglei-full-build\PHengLEIv3d0\Release\PHengLEIv3d0.exe'
# Copy PHengLEI\examples to separate temporary case directories.
# In each copy's bin\cfd_para_hypersonic.hypara, set:
#   maxSimuStep = 2; intervalStepRes = 1;
#   inviscidSchemeName = "roe";        # baseline copy
#   inviscidSchemeName = "entropy_stable"; # new-flux copy
Set-Location <case-directory>
& $exe
```

Both runs exited with code 0 and printed `Job Done.`. The grid contained
79,520 cells. Iteration residuals from the actual solver output:

| Scheme | Iteration | Average residual | Maximum residual |
|---|---:|---:|---:|
| Roe | 1 | `2.16698e-04` | `4.37414e-03` |
| Roe | 2 | `2.16798e-04` | `4.37886e-03` |
| entropy_stable | 1 | `2.16698e-04` | `4.37414e-03` |
| entropy_stable | 2 | `2.09549e-04` | `4.25332e-03` |

These two-step smoke runs verify that both selectors are independently
selectable and complete iterations without a runtime error. The residual
difference after two steps is not evidence of convergence, accuracy, or
stability.

## Unavailable experiments and conclusions

No Sod or Lax case, shock-tube initial-condition option, or corresponding
test configuration was found in the repository's example inputs. The actual
PHengLEI Sod/Lax runs, `L1`/`Linf` errors, convergence orders, oscillation
comparisons, and a full-run discrete-entropy history have therefore **not**
been obtained. The two-step 3dball run also does not provide a discrete
entropy history. No results for these missing experiments are inferred or
fabricated.

The implementation is formula-consistent with the requested division of
responsibility: raw states feed only `f_EC`; Minmod-reconstructed states feed
the entropy-variable dissipation jump and wave speed; Roe remains selectable.
The local same-state dissipation check passes, but the hybrid reconstructed
scheme's global entropy stability remains unproved. The observed sign
mismatches are precisely why it is not described as proven entropy stable.
