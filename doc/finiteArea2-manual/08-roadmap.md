# 8. Roadmap, gaps & maintenance status

[← Developer guide](07-developer-guide.md) · [Index](README.md)

The reduced-order implementation is now the sole `finiteArea` library. This
page tracks remaining functional and maintenance work; it is not a pre-rename
checklist.

## 8.1 Verification gates

| Gate | Status |
|------|--------|
| Regression coverage for geometry and operators | DONE — covered by the current FA tests and tutorial harness. |
| MMS for `grad`/`laplacian`/`div (conv)` | DONE — ~2nd order on refinement. |
| Parallel solver E2E | DONE — `liquidFilmFoam` 4-proc matches serial. |
| Moving-mesh / ALE (`S0`/`S00`/SCL) | DONE (used path) — area machinery + `S0`/`S00` latching to machine ε; `fac::ddt` SCL rate under dilation matches `(1/Δt)(1−1/λ²)` to rel Linf 3.8e-14. Swept-volume `meshPhi` on 2-point edges is degenerate (documented, YAGNI-deferred; no consumer). |
| Performance / memory gate (3-D FV + FA) | Run when changing shared FV geometry or storage; ordinary 3-D FV behaviour remains a required no-regression check. |
| `d2dt2` / `ddt` temporal coverage | Inherited from the FV ddt schemes; add focused FA temporal tests when changing this forwarding path. |

## 8.2 Functional coverage gaps

| Item | Status / plan |
|------|---------------|
| Solvers (`liquidFilmFoam`, `surfactantFoam`, `sphereSurfactantFoam`) | DONE for the supplied tutorial paths. |
| Utilities `makeFaMesh`/`checkFaMesh` | DONE for dictionary-driven generation, persisted-mesh checks, and mesh diagnostics. VTK-specific options remain outside these utilities. |
| Post-processing / conversion (`foamToVTK`, `foamToEnsight`, `src/conversion`, `setFields`) | DONE for the current finite-area paths; retain focused coverage when changing persisted-mesh or surface-writer behaviour. |
| `faMeshSubset` | No separate FA class. Use `fvMeshSubset` on the reduced-order `fvMesh`. |
| Derived BCs (`edgeNormalFixedValue`, `fixedValueOutflow`, `ignore`) | DONE — compiled in `src/finiteArea/Make/files`; standard FV BCs are inherited through the field aliases. |
| Residual scheme consumers (`leastSquaresEdgeInterpolation`) | TODO — still on the removed point-topology API. (`faNVDscheme` is DONE: decoupled from `edgeInterpolation`, which has been deleted; validated via the `liquidFilmFoam` cylinder `Gamma`/`linearUpwind` div schemes.) |

## 8.3 Consumer and link maintenance

All current FA consumers use the single `finiteArea` include directory and
library. Avoid inventing a second `finiteArea2` link or library name in new code.

| Class | Count | Status |
|-------|-------|--------|
| 12.C.1 — redundant dual-linkers (no FA usage, no FA headers) | 8 | DONE — dropped `-lfiniteArea`. |
| 12.C.2a — stale include+link, 0 FA usage | 23 | DONE — dropped the `-lfiniteArea` *link* (stale `-I…/finiteArea/lnInclude` left for the rename). |
| Current FA users | — | Use the single `finiteArea` include directory and link name. |
| 12.C.3 — retire `faDecompose`/`faReconstruct` | — | TODO — blocked only by `reconstructParMesh` (conversion tools now DONE). |
| 12.C.4 — blocked on 12.B / rename | — | Legacy FA solvers/utilities/tests, `reconstructParMesh` FA output, `interfaceTrackingFvMesh` (ALE). (FA-aware `foamToVTK`/`foamToEnsight`/`setFields`/`src/conversion` DONE.) |

## 8.4 Compatibility and public API

| Item | Status |
|------|--------|
| `faSchemes`/`faSolution` dict compat (D2) | DONE (fallback reader); 3-D vs 2-D `system/` dirs confirmed separate. |
| Multi-FA-region / `regionProperties` | Plumbing exists; retain end-to-end coverage for each region and post-processing path. |
| Public API / umbrella headers (`faCFD.H`, `faMesh.H`, `edgeInterpolation`, `faBoundaryMesh`, `faPatch`, …) | Current headers and compatibility forwarders are the supported surface; add coverage when changing aliases or include paths. |
| Tutorials + CI | DONE — all 3 legacy FA tutorials migrated to `tutorials/finiteArea/<solver>/<case>` and wired into `tutorials/finiteArea/Alltest` (run-to-completion + field-parity gates via `test/faFieldDiff.py`, refs in `<case>/validation/`). **7/7 gates PASS.** Tolerances: plane `Cs` legacy-**exact** (`Linf=0`, flat); sphere `Cs` `relL2≤5e-3` (curved-mesh flat-chord vs geodesic `deltaCoeffs`, time-stationary ≈0.2 %) + mass drift `≤1e-6` (obs `0`); film `h` self-regression `relL2≤1e-6`. No case blocked (none need `leastSquaresEdgeInterpolation`/ALE/multi-region). |
| Rename + delete | COMPLETE — the reduced-order implementation is `src/finiteArea` and builds `libfiniteArea`. |

## 8.5 Notable deliberate non-changes

- **`fam::ndiv`** implicit normal projection `n(n&·)` cannot be represented on a
  scalar-coefficient `faMatrix<Type>` (needs a block-coupled matrix); the legacy
  library has the identical passthrough. Left as an explicit passthrough with an
  accurate comment.
- **`edgeInterpolationScheme` kept** — it is not a duplicate of
  `surfaceInterpolationScheme` (embeds the FA tangential transform).
- **`div`/`grad`/`laplacian` scheme consolidation** — measured to be equal to the
  FV operator plus a thin tangential correction, but left as-is: the value is
  modest and the convection-dominated `fam::div` usage is genuinely FA-specific.

## 8.6 Maintenance criteria

Keep the tutorial and focused test coverage passing, preserve ordinary 3-D FV
behaviour, and port residual finite-area consumers where they provide real user
value. In particular, `leastSquaresEdgeInterpolation` still needs a
reduced-order point-addressing port before it can be enabled.
