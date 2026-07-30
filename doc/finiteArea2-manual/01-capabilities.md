# 1. Capabilities

[← Index](README.md) · [Next: Changes from legacy →](02-changes-from-finiteArea.md)

This page summarises what the current `finiteArea` library can do today, and
where implementation work remains. Status labels: **DONE** (implemented + validated),
**PARTIAL** (works on the exercised paths, gaps noted), **TODO** (not yet
available).

## 1.1 Mesh & fields

| Capability | Status | Notes |
|------------|--------|-------|
| Reduced-order 2-D `fvMesh` FA region (`faMesh : public fvMesh`) | DONE | Surface polygons = cells, edges = faces. |
| Area fields (`areaScalarField`, `areaVectorField`, …) | DONE | Aliases of `volField` on the FA mesh. |
| Edge fields (`edgeScalarField`, …) | DONE | Aliases of `surfaceField`. |
| Native reduced-order dimensions | DONE | `S()`≡`V()` (`dimArea`), `Le()` (`dimLength`) — no dimension fix-ups. |
| FA↔parent-volume coupling (`volSurfaceMapping`) | DONE | Coupling addressing owned by `faMesh` (`meshFaces_`, `meshPointLabels_`). |
| Standalone (no parent) area region | PARTIAL | Read/construct works; independent `blockMesh` generation is an additive path (see [roadmap](08-roadmap.md)). |
| Multiple FA regions on one parent | PARTIAL | Plumbing exists (`getAllFaRegionOptions`), not exercised end-to-end. |

## 1.2 Discretisation operators

All of the standard explicit (`fac::`) and implicit (`fam::`) operators are
available and evaluate on the reduced-order mesh through the FV machinery plus an
FA tangential projection. See [Geometry & schemes](04-geometry-and-schemes.md)
for the numerics.

| Operator | Status | Notes |
|----------|--------|-------|
| `fac::grad`, `fam::` gradients | DONE | Gauss, least-squares, limited variants. |
| `fac::div` / `fam::div` (convection) | DONE | `gaussFaConvectionScheme` + NVD/upwind/blended edge schemes. |
| `fac::laplacian` / `fam::laplacian` | DONE | Corrected (non-orthogonal) + uncorrected; **coupled-patch (processor) correction works in parallel**. |
| `fac::ddt` / `fam::ddt`, `d2dt2` | DONE | Forward directly to `fvc`/`fvm` (Euler, backward, CrankNicolson, …) — SCL via `V0/V`. |
| `fac::lnGrad` (surface-normal gradient) | DONE | Corrected/uncorrected/limited/fourth; coupled patches use the `deltaCoeffs` form. |
| `fac::edgeIntegrate`, `facAverage` | DONE | |
| `fam::ndiv` (implicit normal divergence) | PARTIAL | Explicit passthrough (identical to legacy); a true `n(n&·)` implicit term needs a block-coupled matrix. |

Scheme selection is read from the region `system/finite-area/fvSchemes` /
`fvSolution`. Existing `faSchemes`/`faSolution` files remain supported as a
deprecated fallback (see [Migration](06-migration-guide.md)).

## 1.3 Boundary conditions

| BC | Status |
|----|--------|
| All standard FV patch fields (fixedValue, zeroGradient, calculated, coupled, processor, …) | DONE (inherited from FV) |
| `clampedPlate` (FA-specific, now type-generic) | DONE |
| `fixedValueOutflow`, `ignore` (FA-specific, templated) | DONE — genuine ports; validated with `surfactantFoam` (`Cs` outlet/`fixedValueOutflow`, `Cs` inlet/`ignore`), run to exit 0. |
| `edgeNormalFixedValue` (vector) | DONE — thin backward-compat alias over the FV `surfaceNormalFixedValue` (identical on the reduced-order mesh, since `patch().nf()` ≡ the legacy `edgeNormals()`). Prefer `surfaceNormalFixedValue` in new cases. Validated on `Us` inlet via `surfactantFoam`. |

## 1.4 Parallel

| Capability | Status | Notes |
|------------|--------|-------|
| Live parallel FA mesh construction | DONE | Processor patches are genuine `processorPolyPatch`; empty-subdomain deadlock fixed. |
| `decomposePar` (1→N) | DONE | FA fields split with stock `fvFieldDecomposer`. |
| `reconstructPar` (N→1) | DONE | FA fields merged with stock `fvFieldReconstructor`. |
| `redistributePar -reconstruct` / `-decompose` | DONE | Value redistribution via `faMeshProcAddressing` + `parFvFieldDistributor`. |
| `redistributePar` N→M | DONE (two-step) | reconstruct-then-decompose; fused single-invocation is a future convenience. |
| Parallel solver run (E2E) | DONE | `liquidFilmFoam` on 4 ranks matches serial: `h` bit-identical, `Us` to solver tolerance. |

## 1.5 Ported solvers

Located under
[`applications/solvers/finiteArea`](../../applications/solvers/finiteArea)
(the supplied solver targets use their base names):

| Solver | Binary | Status |
|--------|--------|--------|
| `liquidFilmFoam` | `liquidFilmFoam` | DONE — serial + parallel E2E validated. |
| `surfactantFoam` | `surfactantFoam` | DONE — `Cs` transport runs to completion. |
| `sphereSurfactantFoam` | `sphereSurfactantFoam` | DONE — mass conserved exactly on the closed sphere (curved-manifold check). |

## 1.6 Validation coverage

- **Tutorial regression:** the supplied serial and parallel finite-area tutorial
  drivers exercise the current solver paths and compare stored fields.
- **Regression and comparison coverage:** geometry and operators are exercised
  by the finite-area tutorial harness and the FA test applications. Historical
  legacy comparisons remain useful when changing discretisation conventions,
  but there is no legacy FA library to link alongside the current one.
- **Moving-mesh / ALE:** inherited `S`/`S0`/`S00` support follows the FV mesh
  machinery; swept-volume `meshPhi` on 2-point edges remains a documented
  limitation for moving FA consumers.
- **Parallel:** field redistribution exact to machine epsilon; solver E2E matches
  serial.

See [Developer & testing guide](07-developer-guide.md) for how to run these.
