Transient, moving mesh
----------------------
Two turning rotors

cellVolumeWeight interpolation options
-------------------------------------
The volume method uses geometric overlap volumes. Its default weights preserve
constant fields and do not create interpolation extrema. They are not generally
exact for linear fields and do not make the complete overset solver globally
conservative.

Select the method in system/fvSchemes:

oversetInterpolation
{
    method                  cellVolumeWeight;
    weightingScheme         volume;          // default
    donorSelection          legacy;          // default
    donorFallback           legacy;          // default
    minCoverage             0;               // default, no added coverage gate
    minRetainedFraction      0;               // default, no retained-volume gate
    layerDonorMetric         legacy;          // default
    processorMapMethod       AABB;            // default
}

For quality-aware selection, set donorSelection to coverage. This chooses the
largest retained coverage, with the legacy zone priority resolving numerical
ties. minCoverage is the retained overlap divided by the acceptor's geometric
reference volume, not its FV solver volume or its normalized weight sum.
minRetainedFraction is retained overlap divided by the candidate's original
overlap. Both must be finite fractions in [0,1]. Choose nonzero thresholds for
the mesh and check rejected-cell counts; there is no universal safe cutoff.
For example, minCoverage 0.5 can reject weak geometric coverage, but is not a
universal default. Gates apply before fringe admission and after donor removal.
Complete geometric links remain available to hole reachability; coverage failure
alone is not evidence that a cell lies inside a solid. Enabled gates can change
automatic fringe selection and classify failed ordinary receivers as POROUS.

donorFallback legacy retains quality-mode zone retry and previously blanked-cell
recovery. alternateZones also retains/retries candidates under legacy donor
selection. sameZone additionally allows bounded, ready same-zone face-neighbour
IDW for cells requested by the interpolation fringe. Ordinary receivers retain
their selected blend; historical-hole recovery retains full initialization.
These policies never override geometry-cut holes or selected holeLayers/useLayer.
Fallback candidates/promotions cannot be same-update donors for other receivers.
Candidate memory is conditional on quality/retry controls. IDW weights have no
volume units and are never counted as geometric coverage. Recovery can still
fail; MPI-global summaries distinguish retries, ordinary recovery and failures.

layerDonorMetric overlapWeighted removes first-donor ordering from the layer
score using acceptor FV volume divided by the overlap-weighted donor FV volume.
Accumulation uses physical centres/volumes rather than decomposition-dependent
global labels; provisional recovery uses all ready neighbours with IDW weights.
Legacy scoring remains first-donor dependent. Neither mode promises bitwise
MPI identity at numerical layer ties; coincident donor geometry is diagnosed.
The metric is not a directional aspect-ratio error estimate.
processorMapMethod LOD is an
alternative MPI distribution strategy; benchmark it rather than assuming it
is faster than AABB.

For the optional affine correction:

    weightingScheme         linearExact;
    allowNegativeWeights    false;           // default bounded policy
    stencilExpansionRings   1;               // 0 or 1, never an unlimited walk

The correction minimizes the change from volume weights and enforces constant
and first-moment constraints when feasible. It can add eligible face neighbours
from the selected donor zone. Recovery instead starts from IDW weights. A sample
outside the donor-centroid convex hull cannot have both positive weights and
affine exactness. Infeasible or materially negative corrections fall back to
validated volume/IDW weights, with a first warning and MPI-global reason counts.
The per-update receiver summary distinguishes affine-exact and fallback cells.
Even signed weights can be infeasible for the bounded one-ring donor search.
These fallbacks are not affine-exact. The bounded policy rejects a signed
minimum-change solution; it does not run a constrained nonnegative optimization.

allowNegativeWeights true explicitly permits signed corrections. Use it only
when overshoots and loss of field positivity are acceptable. It is not a
recommended default for volume fractions or positive turbulence quantities.
Moment reproduction does not establish second-order PDE convergence or global
conservation. A partial interpolation blend also depends on the calculated
contribution.

Overlap geometry and restart compatibility
------------------------------------------
The clipping engine triangulates each face from its first vertex and joins
those triangles to the FV cell centre. On non-planar faces, its represented
volume G can differ from FV volume V. Pair bounds, small-overlap acceptance and
coverage now use G. Raw pair volumes and corrected moments are not independently
rescaled, and V is never overwritten. Owned-cell G/V mismatches are reported
globally; mapper debug counts include distributed copies and are labelled local.
References are rebuilt on updates/motion, not cached across geometry epochs.

Fan tets must have positive oriented volume; exactly collinear face triangles
are skipped consistently with clipping. Inverted, non-finite or numerically
unresolved geometry is an explicit error. Roundoff allowances depend on scalar
precision, fan count and coordinate/altitude conditioning. Positive excess is
not blanket-clamped. Material pair excess is an error; enabled candidate-quality
checks reject coverage overcount instead of hiding it through normalization.
A warped-face fixture gives G=1.0133333333 versus V=1.01; the former FV-based
bound falsely rejected it. Changing the first face vertex can change warped-face
triangulation, so this is consistency with that representation, not a unique
continuum volume or a new global FV-conservation guarantee.

Finite restart cellTypes values are rounded to the nearest named classification
and clamped in scalar space before label conversion. Global warnings count unique
affected cells and rounding/clamping categories; integer-sized roundoff is quiet.
Non-finite data or missing files on some non-empty ranks are errors. Missing
optional data on empty ranks is supported. Arithmetic mapFields interpolation of
categorical cellTypes is not physical hole reconstruction, despite compatibility
repair permitting the restart to proceed.

Geometry fixes can change formerly incorrect results: contained/disconnected
overlaps are now found, patch priority is independent of donor order, tiny valid
physical volumes normalize safely, and remote donor demotions propagate to all
dependent stencils. Non-interpolated cells retain the existing HOLE/POROUS
failure classifications. Cells in these classes have no interpolation stencil.

The test applications oversetCellVolumeWeight and oversetCellRecovery cover
analytic overlaps, affine moments, interpolation refinement, MPI/renumbering,
and moving recovery. DebugSwitch cellVolumeWeight reports maximum-rank update
and overlap timings. cellVolumeWeightMethod reports local candidate work and
geometry sizes.
Its diagnostics also quantify accepted pairs and discarded overlap volume.
Debug timing also separates reference-volume work, region reachability, donor
policy, affine work, and front/donor map construction/distribution, with halo
sizes and fixed-point rounds. Phase maxima can come from different ranks and
are not additive.
The tetrahedral geometry engine still determines supported cell geometry;
these changes do not establish robustness for arbitrary warped/non-convex cells.
Large aspect ratios also amplify clipping roundoff. The thin-box regression
(aspect ratio 2.5e6) uses a precision/aspect-ratio-scaled allowance; it is not
evidence of uniform 1e-10 relative overlap accuracy for all cell shapes.
Motion is covered by regressions; arbitrary topology-changing and periodic
solver configurations have not been certified by these tests.

Representative post-correctness profiling
----------------------------------------
Controlled 3-D cases covered two-zone ordinary fallback (2240 cells) and
three-zone alternate retry (3571 cells), volume/signed affine and selected
bounded affine policies, and AABB/LOD on 1/2/4 ranks. There were 20 configurations,
three independent runs each and two updates per run. The following are mean
maximum-rank update times in seconds, not comparisons with the old implementation:

    Workload                 Cells  Ranks   AABB       LOD
    ordinary, volume         2240     1     0.130535   0.129761
    ordinary, volume         2240     4     0.036797   0.037153
    alternate, volume        3571     1     0.491316   0.487333
    alternate, volume        3571     4     0.129570   0.134191
    alternate, signed        3571     4     0.131842   0.134597
    ordinary, bounded        2240     4     0.036621   0.037440

Overlap construction accounts for about 95-99% of update time in this matrix.
Reference-volume work, region propagation and instrumented map work are small;
this does not justify a map/region redesign on these cases. AABB/LOD differences
are workload-dependent, not evidence that one is universally faster. Exact
candidate tests were unchanged across the tested distributions: 11664 per
two-update ordinary run and 38448 per alternate run. Donor validation took at
most three rounds for ordinary fallback and one for the alternate cases.
The three-zone four-rank final donor maps used 474 summed halo slots.

For the four-rank signed/AABB case the largest per-process peak RSS was 88552 KiB;
the maximum sum of rank peaks was 347924 KiB. Shared pages are included, and
rank peaks need not occur simultaneously: this is not aggregate physical memory.
These modest synthetic meshes do not establish large-3-D scaling, old-to-new
speedups or a measured benefit from an unimplemented donor-superset cache.

Gated follow-ups (review items 6-9)
---------------------------------
No surface cutter, communication-cache replacement or fringe-only clipping is
enabled by this work. Hole reachability remains authoritative, but is not an
independent solid inside/outside classifier. A separate planar closed-cube
probe reproduces an active background receiver wholly inside the solid in
serial and on four ranks: a coarse background spans both solid-wall and overset
halo interactions. This is independent of the warped-volume mismatch.

A future solid classifier must require explicitly selected closed, orientable
solid patches, fluid/solid-side orientation, MPI assembly, motion invalidation,
and defined handling of open/empty/extruded-2D/symmetry/periodic boundaries.
Do not treat every PATCH/wall as a solid. triSurfaceSearch::calcInside maps
UNKNOWN/on-surface states to false; a robust design must retain those states and
report unresolved classification instead of silently treating it as fluid.

Map reuse is safe only while every required donor remains in its addressing.
Alternative zones, IDW and affine expansion can add donors; an upfront superset
must be compared with rebuild costs and unused-halo memory. Sparse region
propagation must include coupled/processor and non-local overset links, monotone
status updates and globally agreed termination, not processor faces alone.
Fringe-only clipping remains blocked until independent hole classification is
validated, with complete search for newly exposed/moved cells and required
reachability/recovery/retry links. These redesigns require separate approval.
