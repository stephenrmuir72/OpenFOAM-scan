# 6. Migration guide

[← Parallel & I/O](05-parallel-io.md) · [Index](README.md) · [Next: Developer guide →](07-developer-guide.md)

This page covers updating a solver or case written for the earlier stand-alone
`finiteArea` implementation to the current `finiteArea` library. The library
name and link name are unchanged; the implementation behind them is not.

## 6.1 Porting a solver

The mechanical recipe (validated verbatim on `liquidFilmFoam`, `surfactantFoam`,
`sphereSurfactantFoam`):

1. **`Make/options`** — use the current library include path and link:
   ```
   -I$(LIB_SRC)/finiteArea/lnInclude
   -lfiniteArea
   ```
   Link only this single FA library together with the normal FV dependencies.

2. **`Make/files`** — choose the executable name appropriate for the application.
   The supplied finite-area solvers use their base names (`liquidFilmFoam`,
   `surfactantFoam`, and `sphereSurfactantFoam`).

3. **`<solver>.C`** — before `#include "createFaMesh.H"` add:
   ```cpp
   const word aRegionName(faMesh::prefix());   // default "finite-area" region
   ```
   `createFaMesh.H` uses the in-memory generate → `write()` → reload pattern so
   schemes are populated.

4. **Field construction** (`createFaFields.H` etc.) usually needs **no change**:
   `aMesh.thisDb()`, `aMesh.Le()`, `aMesh.areaCentres()`,
   `aMesh.faceAreaNormals()`, `linearEdgeInterpolate`, single-arg
   `volSurfaceMapping(aMesh)`, and the `faScalarMatrix`/`fam::`/`fac::` operators
   all exist unchanged. `aMesh.S()` resolves to `aMesh.V()`, so
   `sum(Cs*aMesh.S())` works.

5. **API renames** exposed by the port (rare):
   - `evaluateCoupled<coupledFaPatch>()` → `evaluateCoupled<coupledFvPatch>()`
       (current area patches are `fvPatch`; `coupledFaPatch` does not exist).

No second-library API is involved. A straightforward port normally only needs
the field-header and patch-class adjustments below.

## 6.2 Porting a case

### Scheme / solution dictionaries

Preferred: put scheme data under the region directory as FV
dictionaries:

```
system/finite-area/fvSchemes
system/finite-area/fvSolution
```

These are structurally identical to the legacy `faSchemes`/`faSolution` (same
top-level tables, including `lnGradSchemes`). **Backward compatibility (D2):** if
the region `fvSchemes`/`fvSolution` are absent, the legacy
`system/finite-area/faSchemes`/`faSolution` are read as fallback content and a
one-line deprecation notice is emitted. So existing cases keep working
unchanged; migrate at leisure by renaming the files.

> The 3-D parent region reads `system/fvSchemes`/`fvSolution`; the FA region
> reads `system/finite-area[/<region>]/…`. They never collide.

### Field file headers

Legacy FA field files declare `class areaScalarField` / `edgeScalarField` in the
`FoamFile` header. The current finiteArea stores them as `vol`/`surface` fields,
so change
the header class:

```
areaScalarField  → volScalarField
areaVectorField  → volVectorField
edgeScalarField  → surfaceScalarField
```

(This is exactly how the `tutorials/finiteArea/liquidFilmFoam` case ships.) The
internal data format is unchanged.

### Mesh generation

- Coupled to a parent (film/shell): generate the area region from
  `system/finite-area/faMeshDefinition` (`polyMeshPatches` + `boundary`
  sub-dictionary). For a boundary-layer film the workflow is
  `blockMesh` → `collapse2DMesh <region>` → `checkMesh -region <region>` →
  `setFields -region <region>`.
- `faMeshDefinition` is read under its FA name, so no rename is needed there.

## 6.3 Compatibility surface

The following legacy names resolve unchanged so most downstream code compiles
against the current finiteArea:

- Types: `faMesh`, `areaScalarField`/…, `edgeScalarField`/…, `faMatrix`,
  `faScalarMatrix`, `faPatch`, `faPatchField`, `faePatchField`,
  `faPatchFieldBase`.
- Operators: the full `fac::`/`fam::` set.
- Accessors: `S()`/`S0()`/`S00()`, `Le()`/`magLe()`, `areaCentres()`,
  `edgeCentres()`, `faceAreaNormals()`, `pointAreaNormals()`, `faceCurvatures()`,
  `faceLabels()`/`meshFaces()`, `whichFace()`.
- Selectors: `faMesh::TryNew(...)`, `hasMeshFiles(...)`, `mesh(...)`,
  `Registry(...)`.
- Legacy `#include` names: besides `faCFD.H`/`faMesh.H`/`areaFields.H`/
  `edgeFields.H`/`faMatrices.H`/`fac.H`/`fam.H`/`volSurfaceMapping.H`,
  the forwarders in `compat/` keep `faMatrix.H`, `faScalarMatrix.H`,
  `faMatricesFwd.H`, `faPatch.H`, `faPatchField(s)(Fwd).H`,
  `faePatchField(s)(Fwd).H` and `faBoundaryMesh.H` resolving (see
  [architecture §3.3](03-architecture.md#back-compatibility-forwarding-headers-compat)).
  Not shimmable (consumer change at rename): per-type `…FaPatchField` BC *class*
  names → `…Fv[s]PatchField`; the deleted `edgeInterpolation` MeshObject;
  `faSchemes`/`faSolution` classes; sliced FA fields; `faMeshSubset`/mappers.

## 6.4 Known migration caveats

- **Boundary conditions**: standard FV patch fields work directly. The
   finite-area library also registers `clampedPlate`, `fixedValueOutflow`,
   `ignore`, and the compatibility `edgeNormalFixedValue` type.
- **Post-processing**: FA-aware `foamToVTK`, `foamToEnsight` and `setFields` use
   the current finiteArea library (link `-lfiniteArea` only) and are validated on an FA
  case. They read the **persisted** reduced-order mesh, so the on-disk fields
  must carry `boundaryField` entries matching the definition patches
  (`faBoundary`) — constraint patches (e.g. `symmetry`) need the matching
  patchField type. `src/conversion`/`libconversion` is FA-agnostic.
- **Subset operations**: use `fvMeshSubset` on the reduced-order `fvMesh`; there
   is no separate FA subset class.
- **Residual scheme consumer**: `leastSquaresEdgeInterpolation` is not compiled
   because it still depends on the old point-topology API. The other listed edge
   schemes remain available.
