# Overset compatibility notes

The exploratory overset fixes affect shared `cellCellStencil` code,
`inverseDistance`, `trackingInverseDistance`, `cellVolumeWeight`, and
`leastSquares`, not just tracking searches. Exact patch-triangle hole cutting,
candidate ranking/retry, receiver-specific fringe filtering, and corrected
front-layer selection can change cell partitions and existing case results
without changing dictionary settings. Compare moving-case fields, continuity,
and cell-type counts before adopting the branch.

## Patch voxel budget

`trackingInverseDistance` accepts this entry in `oversetInterpolation`:

```foam
maxPatchVoxelCount 16777216;
```

The positive count limits **each zone's patch grid on each processor**, both
initial `searchBoxDivisions` and refinement. It is not a total memory limit:
other zones, processors, geometry, and search structures consume additional
storage. The existing 16,777,216 default is unchanged. For example, a
`(300 300 300)` grid needs a budget of at least 27,000,000. Invalid dimensions,
over-budget allocation, arithmetic overflow, or exhausted refinement fail
explicitly rather than allocating without a bound.

Patch intersection uses the same tet-face triangulation as `CELL_TETS`,
including OpenFOAM's point-zero fallback and warnings when a tet base cannot
be found. Structurally invalid faces with fewer than three points remain errors.

## Donor policy and initialization

For the inverse-distance methods, `allowHoleDonors true` permits hole members
in ordinary expanded stencils,
including across consecutive mesh updates. Referenced hole cells use the
existing SPECIAL extrapolation path. It does **not** make holes eligible
primary donors or relax overset receiver fringe restrictions.

Previously blanked receivers require initialization from donors whose field
values were ready before the update. Old-hole donors remain excluded from
every retained/expanded recovery stencil and same-zone fallback regardless
of `allowHoleDonors`. Ordinary interpolation also retains the previous-hole
readiness filter unless hole donors are explicitly allowed. Changing a donor's
classification during an update
does not initialize its snapshot field values. Recovery is not staged
through newly exposed cells; without eligible ready donors, the receiver
remains a hole with an explicit warning.

Deliberate holes selected by `holeLayers`/`useLayer` remain protected from
recovery. Volume-weight same-zone fallback observes receiver-specific fringe
eligibility, including communicated neighbours, while ordinary overlap
weights retain their volume-weight semantics.
ID/TID share post-front restoration; all three methods share recovery donor
eligibility and receiver-aware same-zone neighbour discovery. CVW's application
of overlap weights remains separate to avoid replacing them with distance weights.

## Numerical weights and diagnostics

`weightingScheme` defaults to `inverseDistance`. `leastSquares` and
`linearExact` select the shared affine-weight calculation; the standalone
`leastSquares` method uses that calculation as before. Scheme names are
validated once when constructing the stencil.

Normalization validation uses compensated extended-precision accumulation
and a tolerance based on scalar precision, stencil length, and weight
magnitude, capped to reject gross cancellation errors. Finite signed weights
are not clipped. Generated affine numerators share a compensated
extended-precision normalization denominator, avoiding inconsistent sums in
the Schur subtraction and final weights. General supplied-weight validation
is not relaxed or automatically repaired.

Full three-dimensional moments use compensated accumulation and an error
budget that accounts for local solve/Schur cancellation and rounding of
absolute donor/receiver coordinates, including flat and rotated meshes.
The local solve allowance is capped at `8*sqrt(epsilon)` times the stencil
span. Coordinate uncertainty amplified beyond that local cap plus the
unamplified input-rounding floor is rejected rather than used to justify an
arbitrarily large residual. This preserves off-plane rejection and explicit
fallback for numerically unresolved stencils without changing the SVD cutoff.
Numerical failure injects from the **nearest eligible finite donor**,
not the first stencil entry; this fallback is normalized but not affine-exact.
Exact matches still inject at the matching donor.

Each stencil instance warns on the first nearest-donor fallback and first
invalid-stencil occurrence. Subsequent affected updates print concise Info
summaries; event-free updates add no summary. Nearest-fallback counts are
MPI-global **weight-generation calls**, including intermediate, discarded,
and recomputed stencils, not unique final receivers. Counts reset after
each report. Invalid-stencil **event counts** report empty/misaligned,
nonfinite, and normalization failures separately from remaining demotions and
calculated-cell adjacency. Recovery tries ready-donor fallback where possible.
Fatal errors, unrecoverable-cell warnings, and debug cell identification are
unchanged.

Affine exactness does not imply boundedness or cell-average conservation;
valid extrapolation weights can be negative. No conservative transfer
algorithm or new weight limiter is introduced.

## Validation and deferred performance work

Focused runners under `applications/test` cover patch intersections,
donor fallback, and moving recovery. They include negative tet-base fallback,
budget/scheme dictionary failures, repeated-update hole-donor policy,
poison-independent recovery, signed-weight validation, nearest injection,
volume-weight fallback eligibility, MPI, and random renumbering as applicable.
The patch-intersection runner itself is serial; moving and donor tests provide
parallel integration coverage. Affine checks additionally retain concrete
translated/rotated planar, anisotropic normalization, and near-donor Schur
cancellation counterexamples, test signed extrapolation and off-plane
rejection, and exercise repeated reporting with events on only one MPI rank.
The moving runner tests 20 fixed updates of `twoSimpleRotors` for ID/TID with
`linearExact` and standalone `leastSquares`, plus MPI TID: zero fallback calls,
finite pressure/velocity probes, and matching global cell partitions.

After stencil class-layout changes, rebuild overset solver consumers as well
as `liboverset`. The moving runner rebuilds `overPimpleDyMFoam` before use;
other overset solvers must likewise match the library's headers.

Large genuinely 3-D retry benchmarking is deferred to user testing.
One-cell-thick empty-patch rotor tutorials are not substitutes for it.
Measure matched repeated runs while varying overlapping-zone count at
comparable cell count, and distinguish flood-fill/front/communication costs
from search construction and solver cost. Optimize only if repeatable median
overhead attributable to retries reaches 10% beyond run-to-run variation.
The earlier 23% tutorial difference is not an attribution to retries.

### Matched rotor comparisons

The library built without compiler warnings in `linux64GccDPInt32Opt`.
All three regression runners passed, including actual MPI donor/recovery
execution and random-renumbered donor checks.

Isolated builds of reference `2f53d8246a`, review `af47bc3805`, and these fixes
were compared using the same solver/compiler/precision. Both incompressible
rotor cases ran 40 fixed steps (`deltaT=0.000125`, `endTime=0.005`), with
64 x 64 x 1 patch grids, the tutorial's donor/layer settings, and output at
step 40. These are matched smoke comparisons, not production-case validation.

| Case | Method | Reference holes | Review/fixed holes | Reference cumulative continuity | Review/fixed cumulative continuity |
|---|---|---:|---:|---:|---:|
| simpleRotor | inverseDistance | 506 | 402 | 2.056e-4 | -5.624e-4 |
| simpleRotor | trackingInverseDistance | 497 | 402 | 8.444e-5 | -5.624e-4 |
| simpleRotor | cellVolumeWeight | 486 | 486 | 3.957e-4 | 3.957e-4 |
| twoSimpleRotors | inverseDistance | 853 | 596 | -6.265e-4 | -1.800e-4 |
| twoSimpleRotors | trackingInverseDistance | 845 | 596 | -3.596e-4 | -1.800e-4 |
| twoSimpleRotors | cellVolumeWeight | 753 | 753 | -2.189e-4 | -1.637e-4 |

For default weighting, final calculated/interpolated/hole counts and serial
cumulative continuity exactly matched review HEAD. Final pressure and velocity
probe series in `twoSimpleRotors` also matched review HEAD. All six final cases
completed with three MPI ranks and the same global cell partitions as serial;
parallel continuity values are not bitwise identical to serial. Additional
three-rank review-HEAD runs of all three methods in `twoSimpleRotors` reproduced
the final parallel partitions, continuity, and all 40 pressure/velocity samples
exactly, confirming those serial/parallel numerical differences predate this fix.

Before the residual-validation follow-up, the three affine variants
(ID/TID with `linearExact`, and `leastSquares`)
completed `simpleRotor` with 402 holes and cumulative continuity -2.628e-5.
They reported nearest-donor fallback during every update; this is explicit
loss of affine exactness for rejected stencils, not a conservation guarantee.
These measurements predate the coordinate/cancellation-aware checks above.
The reviewed adaptive 20-update `twoSimpleRotors` reproduction now has zero
fallback calls and the same cell partitions as before the follow-up. The
permanent fixed-update serial/MPI regressions also complete without fallback;
solver completion alone is not the affine-exactness check, which is covered
by the direct nonconstant scalar/vector fixtures.

Single-run serial solver costs were 2.39/2.33/8.66 seconds for
simpleRotor ID/TID/CVW and 6.98/6.93/20.39 seconds for twoSimpleRotors.
Compared with review HEAD, these costs ranged from 0.983 to 1.017 times the
review measurements. Shared-machine variation and lack of repeated runs
preclude a performance conclusion; the deferred 3-D retry benchmark is still
required before considering optimization.

The Keysight Technologies attribution in `patchFaceIntersection.C/.H` was
confirmed by the contributor and is unchanged.
