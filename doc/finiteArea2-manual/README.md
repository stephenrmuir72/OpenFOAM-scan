# finiteArea Developer & Maintainer Manual

The current `finiteArea` library is the reduced-order implementation developed
under the temporary name `finiteArea2`. It has replaced the earlier
stand-alone finite-area implementation under the existing library name.
`finiteArea` is an implementation of the OpenFOAM finite-area (FA) method in
which an FA mesh is expressed as a **reduced-order 2-D `fvMesh`**: the surface
polygons become finite-volume "cells" and their bounding edges become
finite-volume "faces" (a line — a face of reduced topological order). FA thereby
reuses the entire finite-volume (FV) container, addressing, matrix/solver,
interpolation, boundary and parallel infrastructure, and keeps only the
genuinely different 2-D-manifold *numerics*.

This manual is aimed at developers and maintainers working on or against the
current library. The refactor history is available in the commit history; this
directory does not describe a second installable library.

## Contents

| Page | Purpose |
|------|---------|
| [1. Capabilities](01-capabilities.md) | What the current library can do; feature/coverage matrix; supplied solvers. |
| [2. Changes from the earlier `finiteArea`](02-changes-from-finiteArea.md) | Side-by-side of the old stand-alone library and the current reduced-order design. |
| [3. Code architecture](03-architecture.md) | Layering, the `faMesh : public fvMesh` decision, type aliases, key classes. |
| [4. Geometry & discretisation schemes](04-geometry-and-schemes.md) | `faGeometryScheme`, 2-D primitives, `fac`/`fam` operators, `lnGrad`, tangential transforms. |
| [5. Parallel execution & I/O](05-parallel-io.md) | `decomposePar`/`reconstructPar`/`redistributePar`, `faMeshProcAddressing`, coupled patches. |
| [6. Migration guide](06-migration-guide.md) | Updating users of the earlier stand-alone finite-area API and cases. |
| [7. Developer & testing guide](07-developer-guide.md) | Build layout, test apps, validation harnesses, how to add a scheme/BC. |
| [8. Roadmap, gaps & maintenance](08-roadmap.md) | Remaining functional gaps and maintenance work. |

## The idea in one picture

```
  Earlier finiteArea (stand-alone)      Current finiteArea (thin layer on FV)
   ───────────────────────────          ─────────────────────────────
   faMesh (own topology/geometry)   →    faMesh : public fvMesh   (2-D fvMesh)
   areaField / edgeField            →    volField / surfaceField  (aliases)
   faMatrix / faScalarMatrix        →    fvMatrix                 (alias)
   faPatch / faPatchField           →    fvPatch / fvPatchField   (aliases)
   edgeInterpolation (own deltas)   →    surfaceInterpolation weights + FA transform
   fac::/fam:: (own assembly)       →    fac::/fam:: (FV assembly + tangential projection)
   faMeshDistributor                →    redistributePar + fvMeshTools addressing
```

## Key facts to orient a new maintainer

- **An area field *is* a `volScalarField`; an edge field *is* a
  `surfaceScalarField`** on the FA `fvMesh`. This is why stock FV field
  reconstruct/decompose, matrix assembly and linear solvers work unchanged.
- **`faMesh::S()` is `fvMesh::V()`.** On the reduced-order mesh the polygon
  "cell volume" is the finite-area face area (`dimArea`), and the edge "face
  area" `Le()`/`magLe()` is the in-plane edge normal times edge length
  (`dimLength`). These native dimensions make `fvMatrix` dimension checks sound.
- **2-D geometry lives in one place**: the runtime-selected
  `faGeometryScheme` (`: public basicFvGeometryScheme`), auto-selected when
  `primitiveMesh::nTopologicalD() == 2`.
- **FV backward compatibility is a hard constraint.** Every change to the shared
  FV core (`primitiveMesh`, `fvMesh`, geometry schemes) must be behaviourally and
  performance-neutral for ordinary 3-D FV runs.

## Source locations

- Library: [`src/finiteArea`](../../src/finiteArea)
- Shared FV support changes: [`src/finiteVolume`](../../src/finiteVolume),
  [`src/OpenFOAM/meshes/primitiveMesh`](../../src/OpenFOAM/meshes/primitiveMesh)
- Ported solvers:
  [`applications/solvers/finiteArea`](../../applications/solvers/finiteArea)
- Test / validation apps:
  [`applications/test`](../../applications/test) (FA construction, geometry and
  redistribution tests, plus the tutorial regression drivers)

## Current library identity

There is one FA library: `libfiniteArea`, built from
[`src/finiteArea`](../../src/finiteArea). Consumers use
`-I$(LIB_SRC)/finiteArea/lnInclude` and `-lfiniteArea`. There is no
`libfiniteArea2`, `src/finiteArea2`, or second FA runtime-selection table.

What changed for users of the FA API:

- **`faMesh` is a reduced-order `fvMesh`** (`faMesh : public fvMesh`): surface
  polygons are FV "cells", bounding edges are FV "faces".
- **`area*`/`edge*` are FV typedefs, not distinct types**: `areaField ==
  volField`, `edgeField == surfaceField`, `faMatrix<T> == fvMatrix<T>`, and
  `faPatchField<T> == fvPatchField<T>` (`faePatchField<T> == fvsPatchField<T>`).
- **`faSchemes`/`faSolution` fall back to `fvSchemes`/`fvSolution`** when absent,
  so a dedicated FA dictionary is optional (D2 dictionary compatibility).

**Back-compatibility surface.** The header-only
forwarders in [`src/finiteArea/compat`](../../src/finiteArea/compat) keep legacy
`#include "faPatch.H"`, `"faPatchFields.H"`, `"faMatrix.H"`, etc. resolving.
These are a permanent, supported compatibility layer.

**Four source changes that a compatibility shim cannot cover** — third-party or
out-of-tree code that hits any of these must be edited:

1. Concrete `...FaPatchField` derived class names (a bespoke boundary condition
   must derive from the `fvPatchField` hierarchy, not a separate FA one).
2. Code relying on FA-specific *sliced* field types.
3. Direct references to the old dedicated `faSchemes`/`faSolution` C++ classes
   (now the `fvSchemes`/`fvSolution` machinery).
4. `faMeshSubset` / FA field mappers / `MapFaFields` (subset + mapping now go
   through the FV `fvMeshSubset` / mapping infrastructure).

**Known open follow-up:** the
`interfaceTrackingFvMesh` free-surface tracker has a residual bulk-kinematic
numerical instability; it compiles, links and runs on `libfiniteArea` and is
tracked separately from the library rename.
