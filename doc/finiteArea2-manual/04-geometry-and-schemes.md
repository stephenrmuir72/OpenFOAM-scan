# 4. Geometry & discretisation schemes

[← Architecture](03-architecture.md) · [Index](README.md) · [Next: Parallel & I/O →](05-parallel-io.md)

## 4.1 Reduced-order geometry

On the FA `fvMesh` the topological dimension is 2 (`primitiveMesh::nTopologicalD()
== 2`). The mapping between FV and FA quantities is:

| FV quantity | FA meaning | Dimensions |
|-------------|-----------|------------|
| cell `V()` | face area `S()` | `dimArea` (m²) |
| face area vector `Sf()` | edge integration vector `Le()` (in-plane outward edge normal × edge length) | `dimLength` (m) |
| cell centre `C()` | `areaCentres()` | m |
| face centre `Cf()` | `edgeCentres()` | m |
| `deltaCoeffs()` | edge 1/|d| | 1/m |

`Sf()` and `Le()` are the same vector on the reduced-order mesh (`Sf() ≡ Le()`,
`magSf() ≡ magLe()` to machine precision), so the base Green–Gauss assembly for
`grad`/`div` is identical between FV and FA; only the tangential projection and
curved-manifold transforms are FA-specific.

Additional 2-D members exposed on `fvMesh`: `faceAreaNormals()`,
`edgeAreaNormals()`, `pointAreaNormals()`, `faceCurvatures()`.

## 4.2 The geometry strategy: `faGeometryScheme`

2-D geometry is computed by a dedicated runtime-selected strategy
[`faGeometryScheme`](../../src/finiteVolume/fvMesh/fvGeometryScheme/fa/faGeometryScheme.C)
(`: public basicFvGeometryScheme`, `TypeName("fa")`):

- Auto-selected when `nTopologicalD() == 2` (or explicitly via the `geometry`
  entry in `fvSchemes`).
- `movePoints()` calls the base 3-D `movePoints()` first, then adds the
  2-D-manifold geometry (point/face/edge area normals, `Le`, curvature).
- `basicFvGeometryScheme` is therefore **pure 3-D** again — new geometry
  behaviour is added by adding a scheme, not editing an existing one (SOLID
  Open/Closed).

This is the single home of the algorithm that used to live in the legacy
`faMeshDemandDrivenData.C` and (transiently) in a `basicFvGeometryScheme` 2-D
branch.

### Caching

The tangential-transform cache `edgeTransformTensorsPtr_` is reset in
`fvMesh::clearGeomNotOldVol()` so a stale transform cannot survive a geometry
rebuild (moving mesh).

### Parallel point-area-normals (curved/closed manifolds)

`calcPointAreaNormals` accumulates each point's area-normal locally and then must
**sum the contributions of every rank that shares the point**. Doing this through
the reduced-order mesh's *own* point coupling is wrong: the reduced mesh is
edges-as-faces, so its processor-patch/global point addressing is *collapsed* —
a shared point is treated as coupled to every processor-boundary point. On a
closed manifold this corrupts the normal at **multiply-shared** points: a sphere
pole sits at a `2×2×2` domain corner shared by four *non-adjacent* ranks, and a
unit-seed point sync there returns `30` (≈ the local boundary-point count) rather
than the true sharing count `4`. The summed pole normal collapses into the
tangent plane, which propagates through `fDir1`→`edgeAreaNormals`→`Le` and makes
`fac::grad` (and hence the corrected-scheme Laplacian) parallel-inconsistent —
`max|grad(1)|` was ≈1.9 in parallel vs ≈0.01 serial.

The fix routes the sum through the **parent volume mesh's** global point coupling
(`globalMeshData` `globalPointSlaves`), which is correct and resolves
adjacent-rank *and* multiply-shared points in one collective — the parent-global
analogue of the legacy `faMesh::calcPointAreaNormals` processor exchange +
`sharedPoints` `combineReduce`. The seam is deliberately a **registry lookup**, so
that `basicFvGeometryScheme` (in `libfiniteVolume`) stays generic and the `fvMesh`
ABI is untouched:

- `faMesh::registerParentPointMap()` publishes, in the reduced mesh's own
  registry, a `labelIOList` named `pointAreaNormalsParentMap` = the parent-mesh
  point label of each reduced-order point.
- `basicFvGeometryScheme::calcPointAreaNormals` looks that up and, when present
  (with a `parentMesh()`), uses the subset-map
  `syncTools::syncPointList(parentMesh(), parentMap, result, plusEqOp, …)` instead
  of `syncTools::syncPointList(mesh_, …)`. A genuine (non-FA) `fvMesh`, having no
  such entry, keeps the mesh's own sync; in serial the sync is a no-op.

The parent's `globalData` is independent of the reduced mesh and already built, so
the sum is safe from the construction-time geometry build — no re-entrancy into
`processorPolyPatch::initGeometry` and no floating-point exception on the flat
build. Post-fix, parallel `grad(1)` is bit-identical to serial (≈0.00995), the
poles read radial `≈(0,0,±1)`, and the closed-sphere `Gauss linear corrected`
transport matches serial to solver-ordering tolerance (`relL2≈7e-6`). What remains
for a full `interfaceTrackingFvMesh` port is the *moving-mesh* coupling (old-time
area `S0`/`S00`, `meshPhi`, free-surface point motion); the geometry-side point
coupling those depend on is now correct.

## 4.3 Operators = FV assembly + tangential projection

Every `fac::`/`fam::` operator reuses the FV Green–Gauss / LDU assembly and adds
the FA tangential projection `(I − nn)` where the manifold demands it:

- **`grad`** — interior `fac::grad` is bit-identical to `fvc::grad`; only a
  tangential *boundary* correction is FA-specific.
- **`div`** — for a scalar result the `(I − nn)` projection is a no-op, so
  explicit `fac::div(vf)` equals `fvc::div`. The dominant in-tree use is the
  convection form `fam::div(flux, field)`, handled by `gaussFaConvectionScheme`
  (NVD/upwind/blended edge interpolation + tangential projection).
- **`laplacian`** — the implicit part is identical to `fvm::laplacian` (same
  `nonOrthDeltaCoeffs`, `magSf == magLe`); the non-orthogonal correction uses a
  tangential gradient.
- **`ddt`/`d2dt2`** — forward directly to `fvc`/`fvm`. The Space Conservation Law
  is inherited for free through `EulerDdtScheme`'s `Vsc0/Vsc` (`V0/V`).

These equivalences are quantified by the MMS and legacy-parity harnesses (all to
~machine epsilon / 2nd-order — see [testing](07-developer-guide.md)).

## 4.4 Surface-normal gradient (`lnGrad`) and coupled patches

`fac::laplacian`'s non-orthogonal correction calls `fa::lnGradScheme::lnGrad`,
which sets the edge-normal gradient per patch. The important subtlety is
**coupled (processor) patches**:

```cpp
// src/finiteArea/.../lnGradSchemes/lnGradScheme/lnGradScheme.C
forAll(vf.boundaryField(), patchI)
{
    const fvPatchField<Type>& pvf = vf.boundaryField()[patchI];
    if (pvf.coupled())
        ssfb[patchI] = pvf.snGrad(tdeltaCoeffs().boundaryField()[patchI]);
    else
        ssfb[patchI] = pvf.snGrad();
}
```

The no-argument `snGrad()` is a `NotImplemented` stub on `coupledFvPatchField`, so
coupled patches must use the `deltaCoeffs` form (neighbour-exchange
`deltaCoeffs*(neighbour − internal)`), exactly as `fv::snGradScheme` does. This
mirrors FV and is what makes a `corrected` FA laplacian run in parallel.
`deltaCoeffs(vf)` supplies `nonOrthDeltaCoeffs()` whose processor-patch boundary
values are valid, so no extra geometry is needed.

## 4.5 Edge interpolation and the tangential transform

`edgeInterpolationScheme` is **not** a duplicate of `surfaceInterpolationScheme`:
its `interpolate()` embeds the FA tangential transform (`edgeTransformTensors`)
that rotates vector/tensor fields between adjacent, non-coplanar faces on a
curved manifold. Interpolation *weights* and `deltaCoeffs`, however, come from
`surfaceInterpolation` — FA does not recompute geodesic deltas independently.

Available edge schemes: `linear`, `upwind`, `linearUpwind`, `Gamma`, `blended`,
`skewCorrected`. The `faNVDscheme` used by `Gamma`/`blended`/`linearUpwind` was
the last consumer of the `edgeInterpolation` MeshObject; it now takes its base
weights from `surfaceInterpolation::weights()` and computes the geodesic P→N
distance inline, so the `edgeInterpolation` class has been deleted. One scheme
remains on the residual list: `leastSquaresEdgeInterpolation` (still on the
removed point-topology API) — see [roadmap](08-roadmap.md).

## 4.6 Scheme selection

Scheme data is entered exactly as for FV, in the FA region directory:

```
system/finite-area/fvSchemes      # gradSchemes, divSchemes, laplacianSchemes, ...
system/finite-area/fvSolution
```

For backward compatibility, if `fvSchemes`/`fvSolution` are absent the legacy
`system/finite-area/faSchemes`/`faSolution` are read as fallback content (they are
structurally identical); a one-line deprecation notice is emitted when the
fallback is in force. See [Migration](06-migration-guide.md) for details.
