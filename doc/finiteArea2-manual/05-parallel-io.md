# 5. Parallel execution & I/O

[← Geometry & schemes](04-geometry-and-schemes.md) · [Index](README.md) · [Next: Migration →](06-migration-guide.md)

Because an area field *is* a `volScalarField` and an edge field *is* a
`surfaceScalarField` on the FA `fvMesh`, the stock FV parallel machinery does the
heavy lifting. The only FA-specific piece is deriving the FA region's processor
addressing from the parent decomposition.

## 5.1 `faMeshProcAddressing`

[`faMeshProcAddressing`](../../src/finiteArea/faMesh/faMeshProcAddressing/faMeshProcAddressing.H)
emits the FA region's `{cell,face,point,boundary}ProcAddressing` from the parent
`faceProcAddressing`/`pointProcAddressing` plus the coupling `faceLabels`/
`pointLabels` and the reduced-order topology — by **pure lookups, no
communication**. It replaces the core of the legacy `faMeshReconstructor`
(`calcAddressing`) with no serial-mesh rebuild and no collective half-edge
matching.

| Addressing | Derivation |
|-----------|-----------|
| `cellProcAddressing` (area-face = cell) | `faceLabels` through parent `faceProcAddressing`, ranked into serial order. |
| `pointProcAddressing` | FA point → parent point via `faMesh::pointLabels()`, then parent `pointProcAddressing`. |
| `faceProcAddressing` (edge = face) | local edge identified by its global FA point pair, looked up in the serial FA edge→face map; orientation is the signed `±(serialFace+1)` convention. |
| `boundaryProcAddressing` | identity; `-1` for processor patches. |

Because shared points already carry global identity, processor-boundary edges
(internal in serial) resolve with no extra communication. For SPMD use
(`redistributePar`), the helper stores broadcastable serial key-data.

Validation is provided by the finite-area parallel tutorial driver and by
focused tests of `faMeshProcAddressing`; the addressing contract is a complete
cover of the serial mesh with correctly signed processor-split edges.

## 5.2 `decomposePar` (1→N) and `reconstructPar` (N→1)

Both utilities link the single `finiteArea` library (the legacy
`faDecompose`/`faReconstruct` paths are removed). They:

1. Detect the region via `system/finite-area/faMeshDefinition`.
2. Demand-build the complete FA mesh and per-processor FA meshes.
3. Derive addressing via `faMeshProcAddressing`.
4. Split / merge the FA vol & surface fields with the stock
   `fvFieldDecomposer` / `fvFieldReconstructor` — including coupled processor
   boundary fields.

There is no separate FA mesh decomposition on disk — per-processor FA meshes are
demand-rebuilt by solvers and by `reconstructPar`.

The finite-area tutorial parallel driver validates field split and merge on the
supplied cases; oriented edge fields preserve the signed face-addressing
convention through reconstruction.

> **Edge fields must be `oriented`.** The signed-`faceProcAddressing` flip is only
> applied to a surface field that carries `oriented 1`. Real solver fluxes are
> oriented; synthetic test fluxes must set it explicitly.

## 5.3 `redistributePar`

Links `finiteArea` only. It redistributes
FA field **values** through `faMeshProcAddressing` + `fvMeshTools` maps +
`parFvFieldDistributor` — no FA topology is rebuilt.

- **Reconstruct (N→1)** — build serial + per-proc FA mesh, derive FA proc
  addressing, assemble an in-memory reconstruct `mapDistributePolyMesh` with
  `fvMeshTools::createReconstructMap`, gather field values through
  `parFvFieldDistributor`. Validated 4→1 to machine epsilon.
- **Decompose (1→N)** — the inverse: `fvMeshTools::createDecomposeMap` (flip on
  the send side); the master reads complete fields and scatters to every rank
  (all ranks are writers, `isWriteProc=true`). Cut-edge (processor-patch) values
  are populated from the mapped source values (the `parFvFieldDistributor`
  cut-face fix), so a decompose→reconstruct flux round-trip is exact. Validated
  1→4, including an empty-subdomain `simple (4 1 1)` split (proc 3 = 0 cells, no
  deadlock).
- **N→M** — supported as reconstruct-then-decompose (both exact), i.e.
  `redistributePar -reconstruct` then `redistributePar -decompose`. A *fused*
  single-invocation N→M is a documented future convenience.

Robustness: FA blocks are gated on `-no-finite-area`; parent addressing reads are
size-checked (clear fatal instead of an out-of-range crash).

## 5.4 Live parallel construction & the empty-subdomain deadlock

Processor patches on the reduced-order region are genuine `processorPolyPatch`
built from the parent decomposition. The historically defect-prone
point-area-normal synchronisation and empty-subdomain collective ordering are
handled by `faMesh::finaliseParallelBoundary()`, which drives the boundary
topology/geometry collectives explicitly and uniformly on every rank.

A core deadlock (which affected FV as well) was traced to
`primitiveMesh::nTopologicalD()`'s collective `reduce` first-firing from a lazy,
rank-asymmetric `fvMeshDistribute` geometry path when a subdomain was empty. It is
fixed by seeding `nTopologicalD()` uniformly after topology change
(`fvMesh::updateMesh`) and in `finaliseParallelBoundary()`.

## 5.5 Parallel solver E2E

`liquidFilmFoam` on 4 ranks (`mpirun -np 4 liquidFilmFoam -parallel`) runs to
completion; reconstructing and diffing against the *identical* serial case gives:

- `h`: **bit-identical** (max|Δ| = 0)
- `Us`: relLinf ≈ 2.9e-6, relL2 ≈ 2.7e-7 (linear-solver ordering tolerance)

This required the coupled-patch `lnGrad` fix in
[§4.4](04-geometry-and-schemes.md#44-surface-normal-gradient-lngrad-and-coupled-patches).

The **implicit** analogue was needed for the surfactant solvers, whose transport
equation carries an implicit `fam::laplacian(D, Cs)`: on a processor patch the FA
laplacian assembly called the no-arg `gradientInternalCoeffs()` /
`gradientBoundaryCoeffs()`, an intentional `NotImplemented` stub on
`coupledFvPatchField`. `fa::gaussLaplacianScheme::famLaplacian` now guards
`coupled()` patches and uses the `deltaCoeffs` form, mirroring
`fv::gaussLaplacianScheme`. With that fix `surfactantFoam` runs on 4 ranks and
`Cs` reconstructs to `relL2 ≈ 7.3e-6` vs serial.

This serial+parallel parity is now **continuously gated** by the tutorial CI
harness [`tutorials/finiteArea/Alltest-parallel`](../../tutorials/finiteArea/Alltest-parallel)
(invoked by `Alltest`): `planeTransport` (`Cs`, np=4), `liquidFilmFoam/cylinder`
(`h`, np=4) and `sphereSurfactantFoam/sphereTransport` (`Cs`, np=8; closed
manifold) all PASS (**7/7, no SKIP**); see [§7 → the parallel (MPI) leg](07-developer-guide.md#the-parallel-mpi-leg)
for the exact per-case pipelines and tolerances.

> **Standalone FA region parallel recipe.** `decomposePar` splits the FA
> *fields* but does not persist a per-rank FA *mesh*; a solver that read-builds
> the FA mesh (`createFaMesh.H`, `MUST_READ`) therefore needs the per-rank FA
> mesh persisted by `makeFaMesh -parallel`. Use `decomposePar` **and**
> `makeFaMesh -parallel` (see §7). The film case is different: the solver
> demand-builds the FA mesh from the `collapse2DMesh` region, so no `makeFaMesh`
> is used.

### Closed-surface (closed-manifold) parallel construction

A finite-area region that is a **closed** manifold (e.g. the whole `outer` patch
of a sphere, `sphereSurfactantFoam/sphereTransport`) now **constructs and runs in
parallel**. Two coupled problems had to be resolved in
`faMesh::buildRegionMeshFromDefinition`
([`src/finiteArea/faMesh/faMesh.C`](../../src/finiteArea/faMesh/faMesh.C)):

1. **Neighbour-identity `patchMap`.** The reduced-order processor patches are
   cloned from a parent template patch. The old *positional* map cloned arbitrary
   parent patches; on a closed manifold the region has **no** non-processor
   boundary, so the parent's leading non-processor patch shifted every processor
   patch onto a template with the wrong (or no) neighbour. The map now sends each
   reduced-order processor patch to rank *k* onto the parent `processorPolyPatch`
   to rank *k* (guaranteed to exist — the inter-domain edges were derived from
   exactly that parent patch), and non-processor patches onto distinct
   non-processor parent templates. This feeds the inherited processor-patch
   machinery a **symmetric** comm graph.
2. **Skip the degenerate coupled-face reorder.** `makeMesh` is called with
   `syncParallel=false` so the inherited `reorderCoupledFaces` /
   `processorPolyPatch::order` (a 3D polygonal-face geometric matcher) is **not**
   run on the reduced-order **2-point** processor faces (zero area, undefined
   normal), which it would otherwise rotate — reversing an edge relative to its
   owner area-cell and breaking the cell edge loop. The coupled-face
   correspondence is instead supplied **by construction**: every boundary face is
   emitted sorted by a global, decomposition-invariant parent edge-endpoint key,
   so the two sides already list the shared edges in the same order with opposite
   (outward-from-own-cell) orientation — the processor convention. The boundary
   is finalised uniformly on every rank by `finaliseParallelBoundary()` (§5.4).

The sphere case builds and runs on 8 ranks, **conserves mass**
(`sum(Cs·V)` `relDrift ≈ 8e-6`), and now matches serial in the transported field
(`Cs` `relL2 ≈ 7e-6`); the CI harness gates run-to-completion, mass, **and
field-parity** as **PASS**.

> **Resolved — parallel point-area-normal coupling on a curved manifold.**
> With construction fixed, the closed-sphere `Gauss linear corrected` transport
> was initially decomposition-dependent (reconstructed `Cs` `relL2 ≈ 4e-2` vs
> `≈5e-4` under `uncorrected`), because the non-orthogonal **correction** term is
> the only place the Laplacian consumes `fac::grad(Cs)`, and **`fac::grad` itself
> was parallel-inconsistent**: `grad` of a constant, which must vanish, gave
> `max|grad(1)| ≈ 1.9` (vs `≈0.01` serial). `calcPointAreaNormals` summed each
> shared point's normal through the **reduced-order mesh's own** point coupling,
> which is *collapsed* on the edges-as-faces mesh — a unit-seed point sync at a
> sphere pole returned `30` (≈ the local processor-boundary-point count) instead
> of the true sharing-rank count `4` (the pole is a `2×2×2` domain corner shared
> by four *non-adjacent* ranks), so the pole normal collapsed into the tangent
> plane and corrupted `fDir1`/`Le`/`grad` there. It was **not** a
> `correctedLnGrad`/`TN` correction-path bug — the coupled-edge
> `edgeTransformTensors` `TN`, `nonOrthCorrectionVectors`, `Le` and
> `edgeAreaNormals` are all parallel-consistent.
>
> **Fix.** Route the point-normal sum through the **parent** volume mesh's global
> point coupling (`globalMeshData` `globalPointSlaves`), which resolves
> adjacent-rank *and* multiply-shared points in one collective — the
> parent-global analogue of the legacy `faMesh::calcPointAreaNormals` processor
> exchange + `sharedPoints` `combineReduce`. `faMesh::registerParentPointMap()`
> publishes the reduced→parent point map as a registry `labelIOList`
> (`pointAreaNormalsParentMap`), and `basicFvGeometryScheme::calcPointAreaNormals`
> looks it up and uses the subset-map `syncTools::syncPointList(parentMesh(),
> parentMap, …)` in place of the mesh's own sync (falling back to the own sync for
> a genuine `fvMesh`; a no-op in serial). The registry seam keeps the `fvMesh` ABI
> untouched (see [§4.2](04-geometry-and-schemes.md#parallel-point-area-normals-curvedclosed-manifolds)).
> The parent's `globalData` is independent and pre-built, so the sum is safe from
> the construction-time geometry build (no `initGeometry` re-entrancy / FPE on the
> flat build). Post-fix: parallel `grad(1)` bit-identical to serial (`≈0.00995`),
> poles radial, sphere field-parity `relL2 ≈ 7e-6` (**PASS**). Flat/open regions
> (plane, film) are unaffected (their non-FA path keeps the mesh's own sync).
> What remains for `interfaceTrackingFvMesh` is the moving-mesh coupling
> (old-time area `S0`/`S00`, `meshPhi`, free-surface point motion), no longer
> blocked by the point-normal sum.

## 5.6 Directory layout on disk

FA data lives under the `finite-area` sub-registry, so paths never collide with
the parent 3-D region:

```
<case>/system/finite-area[/<region>]/{fvSchemes,fvSolution,faMeshDefinition}
<case>/constant/finite-area[/<region>]/faMesh/{faceLabels,pointLabels}
<case>/<time>/finite-area[/<region>]/{h,Us,...}
<case>/processorN/<time>/finite-area[/<region>]/...
```

## 5.7 Multiple FA regions in parallel

Several finite-area regions can share one parent decomposition. Each region is
an independent reduced-order `fvMesh` registered under its own name in the single
`finite-area` sub-registry, so the per-`<region>` path split above applies to
every region simultaneously and each rank writes
`processorN/<time>/finite-area/<region>/…`. Nothing FA-specific is needed for the
coupling: each region grows genuine `processorPolyPatch`es
(`procBoundary0to1`, …) across the parent cut and exchanges halos through the
stock `processorFvPatchField`, exactly as the single-region path in §5.4.

The multi-region path should be validated with two regions such as `filmA` and
`filmB` on a 2-way-decomposed parent; both must construct with
processor patches, resolve their own `system/finite-area/<region>/fvSchemes`
(distinct `Gauss linear` vs `Gauss upwind`, neither picking up the parent
sentinel), evaluate operators independently, and round-trip fields to their
per-region `processorN` directories — all assertions PASS in both serial and
2-rank parallel. Reconstruct/redistribute across *multiple* FA regions is a
residual (single-region N↔1 is covered in §5.2/§5.3).

## 5.8 Parallel post-processing output (ensight / VTK)

FA area/edge output — the [`areaWrite`](../../src/functionObjects/utilities/areaWrite)
function object and the `foamToEnsight` / `foamToVTK` FA converters — now runs in
parallel. There is **no FA-specific gather**: because `faMesh::patch()` yields a
`uindirectPrimitivePatch` over the parent mesh, the surface and its fields ride
the **stock `surfaceWriter` parallel merge** (`PatchTools::gatherAndMerge` for the
geometry + `globalIndex::gather` for each field), i.e. the same gather-to-master
the FV ensight/VTK writers use. Ranks whose decomposition holds no FA faces
contribute an empty patch and participate in the collective uniformly by
construction, so there is no empty-subdomain special case.

The enabling fix was **not** in the writer. In the solver the single FA region is
named `faMesh::prefix()` (`"finite-area"`), and
`faMesh::readCouplingAddressing()` / `writeObject()`
([`src/finiteArea/faMesh/faMesh.C`](../../src/finiteArea/faMesh/faMesh.C))
previously mis-treated that name as a *named sub-region*, building a **doubled**
path `constant/finite-area/finite-area/faMesh/faceLabels` instead of the
single-level `constant/finite-area/faMesh/faceLabels` written by
`makeFaMesh -parallel`. The coupling `faceLabels` then never loaded, `meshFaces_`
stayed empty, and `patch()` returned a **0-face** patch — so the surface merged to
zero faces while the field gather still carried the full per-rank field, a
non-uniform-size collective that desynchronised and crashed (`FOAM FATAL IO ERROR
… incorrect first token … word 'd'`). Both functions now treat
`areaName_ == faMesh::prefix()` (alongside `polyMesh::defaultRegion` / empty) as
the **default** region, so the read and write paths agree on `<prefix>/faMesh`.

**Validated** (np=4, `planeTransport`, gated by `Alltest-parallel`): the parallel
`areaWrite` ensight case opens and its geometry is byte-identical to the serial
run; the reconstructed `Cs@6` matches the serial reference at `relL2 ≈ 7.3e-6`
(`Linf ≈ 3.3e-5`); the VTK `.vtp` output is equivalent (same point/polygon
counts). Serial `areaWrite` / `foamToEnsight` are unchanged (serial reconstruct
bit-identical). The planeTransport parallel leg consequently drops the former
`-noFunctionObjects` workaround; the film/sphere legs retain it only because they
ship no FA function object.
