# 2. Changes from the earlier `finiteArea`

[← Capabilities](01-capabilities.md) · [Index](README.md) · [Next: Architecture →](03-architecture.md)

The earlier `src/finiteArea` library (~71k lines) was a nearly self-contained
re-implementation of the FV method for surfaces: it had its own mesh, fields,
matrix, interpolation, boundary and parallel classes that paralleled — and
largely duplicated — the FV equivalents. The current `finiteArea` removes that
duplication by making the FA mesh a specialisation of `fvMesh`.

## 2.1 What changed conceptually

| Concern | Earlier stand-alone `finiteArea` | Current `finiteArea` |
|---------|--------------------|---------------|
| Mesh | `faMesh` — bespoke class with its own topology, addressing and geometry (`faMeshDemandDrivenData.C`). | `faMesh : public fvMesh` — a reduced-order 2-D `fvMesh`; edges are FV faces, polygons are FV cells. |
| Fields | `areaField`/`edgeField` families (own `GeometricField` instantiations). | `typedef volMesh areaMesh; typedef surfaceMesh edgeMesh;` — area/edge fields **are** vol/surface fields. |
| Matrix | `faMatrix`, `faScalarMatrix`. | `using faMatrix = fvMatrix; using faScalarMatrix = fvScalarMatrix;`. |
| Patches | `faPatch`, `faPatchField`, `faePatchField`, `faBoundaryMesh`. | `faPatch`≡`fvPatch`, `faPatchField`≡`fvPatchField`, `faePatchField`≡`fvsPatchField` (aliases). |
| Interpolation deltas | `edgeInterpolation` computed its own geodesic weights/deltas. | Weights/`deltaCoeffs` come from `surfaceInterpolation`; only the FA **tangential transform** is FA-specific. |
| Geometry | Inside the mesh class. | Runtime-selected `faGeometryScheme` strategy (see [geometry](04-geometry-and-schemes.md)). |
| Parallel redistribution | `faMeshDistributor`, `faFieldReconstructor`, `faReconstruct`/`faDecompose`. | Stock `redistributePar`/`decomposePar`/`reconstructPar` + `faMeshProcAddressing` + `fvMeshTools` addressing. |
| Scheme dictionaries | `system/finite-area/faSchemes`, `faSolution`. | `system/finite-area/fvSchemes`, `fvSolution` (with legacy `fa*` fallback — D2). |

## 2.2 What stayed FA-specific (kept, not duplication)

Not everything in the old library was redundant. The following are genuinely
different for a 2-D manifold and are retained in the current `finiteArea`:

- **`fac`/`fam` operators.** The public API is preserved, but each operator is
  now *FV assembly + an FA tangential projection* `(I − nn)`, rather than a fully
  independent implementation.
- **`edgeInterpolationScheme`.** Reclassified as **keep, not duplication**: its
  `interpolate()` embeds the FA tangential transform (`edgeTransformTensors`)
  that rotates vectors/tensors between adjacent non-coplanar faces — something
  `surfaceInterpolationScheme` does not do. (The lighter `edgeInterpolation`
  *MeshObject* has been deleted entirely — its last consumer, `faNVDscheme`, now
  takes weights from `surfaceInterpolation` and computes the geodesic P→N
  distance inline.)
- **`lnGrad` schemes** (`corrected`/`limited`/`fourth`/`uncorrected`): the
  surface-normal gradient with the FA non-orthogonal correction.
- **`volSurfaceMapping`.** The FA↔3-D-volume bridge; needs the parent-face per FA
  cell, which is why `faMesh` is a subclass rather than a pure `typedef` (see
  decision D1 in [architecture](03-architecture.md)).
- **FA-specific BCs** (`clampedPlate`, `fixedValueOutflow`, `ignore`, and the
  compatibility `edgeNormalFixedValue` type).

## 2.3 Removed / replaced subsystems

- **`faMeshDistributor` / `faFieldReconstructor`** → replaced by
  `faMeshProcAddressing` (derives FA `{cell,face,point,boundary}ProcAddressing`
  from the parent decomposition by pure lookups) feeding stock FV
  field (de)composers and `redistributePar`.
- **`genericFaPatchField` / `genericFaePatchField`** → removed;
  `genericFvPatchField` already registers the `generic` type for FA fields, so
  `libgenericPatchFields` no longer links the legacy FA library (this removed an
  ODR/runtime-selection clash that corrupted the heap when a utility linked both
  FA libraries).
- **`edgeFields.C`** → deleted (was a trivial/empty translation unit under the
  reduced-order design).

## 2.4 Shared-core changes (FV-visible, behaviour-neutral)

Because FA now rides on the FV core, a handful of shared classes gained 2-D
awareness. These are constrained to be behaviourally and performance-neutral for
ordinary 3-D FV:

- `primitiveMesh` gained `nTopologicalD()` and 2-D primitives (`cellDir0`,
  `faceDir1`, `faceVec2`, `pointDir0`). `nTopologicalD()` is a **collective
  reduce**; it is seeded uniformly after topology change so it cannot first-fire
  from a rank-asymmetric lazy geometry path (the empty-subdomain deadlock — which
  affected FV `fvMeshDistribute` too — was traced to and fixed here).
- `fvMesh` exposes FA-named accessors over FV storage (`S()`→`V()`, `S0()`→`V0()`,
  `areaCentres()`→`C()`, `edgeCentres()`→`Cf()`) plus real 2-D members (`Le()`,
  `magLe()`, `faceAreaNormals()`, `edgeAreaNormals()`, `pointAreaNormals()`,
  `faceCurvatures()`).
- `fvMeshGeometry` selects native reduced-order dimensions
  (`nTopologicalD()==2 ? dimLength : dimArea` for `Sf`,
  `nTopologicalD()==2 ? dimArea : dimVolume` for `V`).
- `EulerDdtScheme` / `fvmSup` now declare `fvMatrix` dimensions from
  `mesh().V().dimensions()` instead of a hard-coded `dimVol`, so implicit source
  terms are dimensionally correct on the reduced-order mesh (`dimArea`). This is
  a no-op in 3-D where `V().dimensions() == dimVol`.

A performance/memory no-regression check for these shared-core changes remains
part of the maintenance criteria (see [roadmap](08-roadmap.md)).
