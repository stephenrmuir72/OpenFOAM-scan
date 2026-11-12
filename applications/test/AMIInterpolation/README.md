# cyclicAMI / faceAreaWeightAMI regression tests

This directory contains the AMI regression driver, analytical fixtures, the
`solidFoam`/`pureZoneMixture` fixture template, and the standard-library
runner. The tests are separate from production AMI code. S02/S03 record
baseline behavior; S04 validates the post-fix C0 build.

## Run the C0 S02 gate

Run from a clean shell after sourcing the current checkout:

```sh
source /home/pc/gitRepos/OpenFOAM/openfoam/etc/bashrc
export AMI_IMPROVEMENTS_ROOT=/home/pc/OpenFOAM/OF2312/tests/my_work/03_KTM/14_AMI_improvements
applications/test/AMIInterpolation/Allrun
```

`Allrun` builds `Test-AMIInterpolation` in a fresh external build directory,
runs the targeted serial and MPI tests, validates the solid fixture inputs and
meshes, and writes a JSON report plus an input-hash manifest. The default
artifact paths are under `06_build/current-before/S02`,
`04_runs/current-before/S02`, and `05_results/S02` below
`AMI_IMPROVEMENTS_ROOT`. The runner refuses to put generated artifacts inside
the source checkout. Each invocation uses a unique run ID and does not replace
older results.

The runner accepts `runTests.py` options through `Allrun`. For example,
`--timeout 120` changes the per-command timeout. `--tests` can select a
diagnostic subset, but a subset does not satisfy the full S02 exit gate.
`python3 runTests.py --list` prints the runnable test IDs and expected assertion
IDs.

The default S02 mode accepts only the current C0 checkout. The runner checks
the `WM_PROJECT_DIR`, version, build options, tool paths, linked OpenFOAM
libraries, and source-tree identity.

## S03 B0/C0 observation mode

S03 can build the same test driver separately against B0 or C0. Source the
selected target's environment, then invoke `runTests.py` with `--stage S03`,
`--target-id B0` or `C0`, and `--observation-mode`. Use a distinct build, run,
and result root for each target. For example, with `AMI_ROOT` set to the AMI
improvements workspace:

```sh
source /usr/lib/openfoam/openfoam2606/etc/bashrc
python3 applications/test/AMIInterpolation/runTests.py \
  --stage S03 --target-id B0 --observation-mode --build \
  --project-dir /usr/lib/openfoam/openfoam2606 \
  --build-root "$AMI_ROOT/06_build/baseline-v2606/S03" \
  --run-root "$AMI_ROOT/04_runs/baseline-v2606/S03" \
  --results-dir "$AMI_ROOT/05_results/S03" \
  --validate-solid-fixtures --measure-performance \
  --performance-iterations 100 --timeout 180
```

Repeat in a fresh shell for C0, using its own artifact roots:

```sh
source /home/pc/gitRepos/OpenFOAM/openfoam/etc/bashrc
python3 applications/test/AMIInterpolation/runTests.py \
  --stage S03 --target-id C0 --observation-mode --build \
  --project-dir /home/pc/gitRepos/OpenFOAM/openfoam \
  --build-root "$AMI_ROOT/06_build/current-before/S03" \
  --run-root "$AMI_ROOT/04_runs/current-before/S03" \
  --results-dir "$AMI_ROOT/05_results/S03" \
  --validate-solid-fixtures --measure-performance \
  --performance-iterations 100 --timeout 180
```

Observation mode records assertion outcomes as `observed_pass` or
`observed_failure`. It does not apply the S02 rule that known C0 failures are
expected, and it does not treat a failed assertion as a passing regression.
Only the selected AMI test subset and fixture-input checks are covered by this
runner status. It does not by itself complete the S03 stage gate, which also
requires the tutorial, `solidFoam` runtime, and performance evidence.
With `--measure-performance`, it records build and command wall time, GNU
`time` peak RSS, repeated T01 constructor/mapping/interpolation timings, T02
cache-restore/interpolation timings, and the T03 topology rebuild time. MPI RSS
is the launcher process peak, not an aggregate of rank memory.

## S03 application screens

Run the application screens in a separate shell after sourcing the selected
target. Use distinct external run roots for B0 and C0:

```sh
source /usr/lib/openfoam/openfoam2606/etc/bashrc
python3 applications/test/AMIInterpolation/runS03ApplicationScreens.py \
  --target-id B0 --project-dir /usr/lib/openfoam/openfoam2606 \
  --run-root "$AMI_ROOT/04_runs/baseline-v2606/S03/applications" \
  --results-dir "$AMI_ROOT/05_results/S03/applications" \
  --components solid moving-solid mixer
```

For C0, source the current checkout and use its project directory plus
`current-before` run root. The `solid` component runs T14-T17, a conformal T15
control, and T19 in serial and on two and four ranks with GAMG. `moving-solid`
runs T18 with face creation and caching off, then runs point motion without
face creation both with and without caching. The cache-enabled case records
the `AMICache` interpolation indices and fraction. `mixer` runs copied
`mixerVesselAMI2D` and `mixerVesselAMI2D-topologyChange` cases.

The runner records runtime `pureZoneMixture` selection, material properties,
thermal metrics, cell-zone membership, rank count, AMI rebuilds and created
faces, cache reuse and fraction, reconstruction status, warnings, and
source-file hashes. It compares T19 serial fields with the two- and four-rank
fields and compares the T18 cached and uncached point-motion fields.

Mixer cases run only from their copied directories. The runner sets
`endTime 0.002` and `writeInterval 0.001` so the short run writes fields for
reconstruction; it records the tutorial's actual `numberOfSubdomains`. This is
a setup/motion/transfer screen, not a periodic-accuracy acceptance run.

The completed S03 AMI reports are
`$AMI_ROOT/05_results/S03/S03-B0-20261020T090604Z-55ab99d6.json` and
`$AMI_ROOT/05_results/S03/S03-C0-20261020T090640Z-a3906c55.json`. The latest
solid and moving-solid reports are
`$AMI_ROOT/05_results/S03/applications-final/S03-applications-B0-20261020T092440Z-7ded22c8.json`
and
`$AMI_ROOT/05_results/S03/applications-final/S03-applications-C0-20261020T092458Z-01ae9c7a.json`.
The mixer results are in the B0 and C0 reports under
`$AMI_ROOT/05_results/S03/application-screens-final-2/`.

These are baseline observations, not acceptance claims. Both targets restored
the T18 exact cache hit with recorded weight `0`; cached and uncached fields
matched across all 12 cells. The topology-changing T18 run created four
additional faces and preserved the saved cell-zone membership. Its uniform
temperature error was `0.002775057` on B0 and `0.231249996` on C0, while the
point-motion control and cache run were about `8.6e-7`. T19's normalized
maximum serial-to-MPI temperature difference was about `0.13304` on both
targets. The T09 distributed-to-local state transition passed on C0 but failed
on B0 in all four multi-rank scenarios, with stale maps and source/target
coverage errors of `0.234375` and `0.34375`. The T02 test separately exercises
fractional cache interpolation and records the pre-change defect. Preserve
these observations for later implementation and acceptance steps.

## S04 post-fix fractional-cache gate

S04 verifies the `AMICache::weight()` scalar accessor and its consumers. It
builds `libmeshTools`, `libfiniteVolume`, and `solidFoam` outside the C0 source
tree. Do not use `runTests.py --build` as a substitute for the library builds;
that option builds only the test driver.

Run from the C0 repository root in a fresh shell after sourcing C0:

```sh
source /home/pc/gitRepos/OpenFOAM/openfoam/etc/bashrc
set -e
AMI_ROOT=/home/pc/OpenFOAM/OF2312/tests/my_work/03_KTM/14_AMI_improvements
S04_ROOT="$AMI_ROOT/06_build/current-fixed/S04"
PRODUCTS="$S04_ROOT/products"
C0_LIBBIN="$FOAM_LIBBIN"
mkdir -p "$PRODUCTS/build-root" "$PRODUCTS/lib" "$PRODUCTS/bin"
export FOAM_LIBBIN="$PRODUCTS/lib"
export FOAM_APPBIN="$PRODUCTS/bin"
export FOAM_BUILDROOT="$PRODUCTS/build-root"
export LIBRARY_PATH="$FOAM_LIBBIN:$C0_LIBBIN:$C0_LIBBIN/sys-openmpi"
export LD_LIBRARY_PATH="$FOAM_LIBBIN:$C0_LIBBIN:$C0_LIBBIN/sys-openmpi:$LD_LIBRARY_PATH"

(cd "$WM_PROJECT_DIR/src/meshTools" && \
  wmake -build-root="$FOAM_BUILDROOT" libso)
(cd "$WM_PROJECT_DIR/src/finiteVolume" && \
  wmake -j4 -build-root="$FOAM_BUILDROOT" libso)
(cd "$WM_PROJECT_DIR/applications/solvers/heatTransfer/solidFoam" && \
  wmake -j4 -build-root="$FOAM_BUILDROOT")

python3 applications/test/AMIInterpolation/runTests.py \
  --stage S04 --target-id C0 --project-dir "$WM_PROJECT_DIR" \
  --build-root "$S04_ROOT" \
  --run-root "$AMI_ROOT/04_runs/current-fixed/S04" \
  --results-dir "$AMI_ROOT/05_results/S04" \
  --library-root "$PRODUCTS/lib" \
  --solid-foam "$PRODUCTS/bin/solidFoam" \
  --build --timeout 180
```

The runner checks C0 provenance, builds the test driver in a fresh external
directory, and requires T01, all T02 serial/MPI fraction scenarios, all T08
and T09 scenarios, the cache-plus-face-creation rejection, T13, and the T18
`solidFoam`/`pureZoneMixture` field/GAMG screen. The T18 report must record a
fraction with `0 < w < 1` and a GAMG solve at the same time. Available field
metrics prove that the solver produced output; they do not replace the S13
thermal-accuracy acceptance tests.

## S05 post-fix rotational area-frame gate

S05 verifies that topology-created cyclicAMI faces store the neighbour area
in the neighbour frame. It rebuilds only `libmeshTools` outside the C0 source
tree, then runs the T03 rotation variants, T08, and T09 with pass-only
expectations. Serial and two-rank T03 variants each write a `t03Variant`
file, and the runner copies it into every `processor*/system` directory.

Run from the C0 repository root in a fresh shell after sourcing C0:

```sh
source /home/pc/gitRepos/OpenFOAM/openfoam/etc/bashrc
set -euo pipefail
AMI_ROOT=/home/pc/OpenFOAM/OF2312/tests/my_work/03_KTM/14_AMI_improvements
S05_ROOT="$AMI_ROOT/06_build/current-fixed/S05-final"
PRODUCTS="$S05_ROOT/products"
C0_LIBBIN="$FOAM_LIBBIN"
mkdir -p "$PRODUCTS/build-root" "$PRODUCTS/lib"
export FOAM_LIBBIN="$PRODUCTS/lib"
export FOAM_BUILDROOT="$PRODUCTS/build-root"
export LIBRARY_PATH="$FOAM_LIBBIN:$C0_LIBBIN:$C0_LIBBIN/sys-openmpi"
export LD_LIBRARY_PATH="$FOAM_LIBBIN:$C0_LIBBIN:$C0_LIBBIN/sys-openmpi:${LD_LIBRARY_PATH:-}"

for lib in "$C0_LIBBIN"/*.so; do
  name="$(basename "$lib")"
  if [ "$name" != "libmeshTools.so" ] && [ ! -e "$FOAM_LIBBIN/$name" ]; then
    ln -s "$lib" "$FOAM_LIBBIN/$name"
  fi
done

(cd "$WM_PROJECT_DIR/src/meshTools" && \
  wmake -build-root="$FOAM_BUILDROOT" libso) \
  2>&1 | tee "$S05_ROOT/rebuild-meshTools.log"

python3 applications/test/AMIInterpolation/runTests.py \
  --stage S05 --target-id C0 --project-dir "$WM_PROJECT_DIR" \
  --build-root "$S05_ROOT" \
  --run-root "$AMI_ROOT/04_runs/current-fixed/S05-final" \
  --results-dir "$AMI_ROOT/05_results/S05" \
  --library-root "$FOAM_LIBBIN" --build --timeout 600
```

The runner requires `FOAM_LIBBIN` to equal `--library-root`, and the library
root must lie under `--build-root`. It checks that the driver resolves
`libmeshTools.so` to that root. The S05 stage has no expected failures. Any
failure, timeout, or missing scenario makes the stage fail.

The S05 gate requires 32 passing scenarios: 13 T03 variants (11 serial
and 2 with two ranks: 90 degrees and 53.13 degrees about z), 13 T08 topology
scenarios, and 6 T09 state-transition scenarios including the
cache-plus-face-creation rejection. The serial T03 matrix includes 0, ±30,
90, and 180 degree rotations, an additional ±53.13 degree pair, an arbitrary
axis and non-zero centre, a zero-separation translational control, and a
non-zero translation with `separationVector (0 2 0)`. The input manifest
records the production source and rebuilt `libmeshTools.so` hashes.

Do not run S02 in C0 expected-failure mode against an S05 library. The fix
removes the T03 defect, so S02 is expected to report `T03_NEIGHBOUR_AREA_FRAME`
as an unexpected pass. This is inferred from the code and was not run. S02
remains the pre-fix record. The runner refuses S02 with
an external library root, so this combination is blocked.

## S06 post-fix low-weight frame gate

S06 verifies that low-weight default values stay in the receiver frame on
rotational cyclicAMI faces. It builds the test application from the current
source and runs T04 (serial and two-rank distributed, with pass-only
expectations), T08, and T09 with the cache-plus-face-creation rejection. The
two-rank T04 case is decomposed by the runner with the rotation fixture's
`system/decomposeParDict` (two subdomains, `simple`, `(2 1 1)`), and only
`controlDict` is copied into each `processor*/system` directory.

Run from the C0 repository root in a fresh shell after sourcing C0:

```sh
source /home/pc/gitRepos/OpenFOAM/openfoam/etc/bashrc
AMI_ROOT=/home/pc/OpenFOAM/OF2312/tests/my_work/03_KTM/14_AMI_improvements
python3 applications/test/AMIInterpolation/runTests.py \
  --stage S06 --target-id C0 --project-dir "$WM_PROJECT_DIR" \
  --build-root "$AMI_ROOT/06_build/current-fixed/S06-final" \
  --run-root "$AMI_ROOT/04_runs/current-fixed/S06-final" \
  --results-dir "$AMI_ROOT/05_results/S06" --build --timeout 900
```

The S06 gate requires 21 passing scenarios: 2 T04 variants (serial and two
ranks), 13 T08 topology scenarios, 5 T09 state-transition scenarios, and the
cache-plus-face-creation rejection. The four-rank T04 run is not in the
runner. It was run by hand with a four-subdomain `decomposeParDict` and
passed; the evidence is in `05_results/S06-T04-extension/`. The S06 stage
has no expected failures.

## Result interpretation

The C0 defect tests must fail their named assertions. The runner records these
as `expected_failure` only when the exact declared assertion IDs fail and the
process exits nonzero. An unrelated failure, timeout, missing assertion,
unexpected pass, non-finite value, or changed assertion set is a suite
failure. `expected_failure` means that the regression detected the known C0
defect; it does not mean the defect is acceptable.

| ID | S02 execution |
|---|---|
| T00 | Timeout process-group self-check and solid-reference unit tests |
| T01 | Planar nonmatching transfer control and weighted-adjoint check |
| T02 | Cached midpoint fraction; expects the current C0 assertion failures |
| T03 | 90-degree topology-created face area; expects the current C0 frame failure |
| T04 | Rotational low-weight fallback; checks every face for U (vector), S (scalar control), T (tensor) and ST (symmTensor): receiver-side defaults stay unrotated, donor-valued faces match each side's `forwardT` (S06 fix; passes on the fixed build, fails pre-fix for U, T and ST) |
| T05 | Seed selection with and without restart: every face ordering (`as-built`, `src-reversed`, `tgt-reversed`, `both-reversed`) must give the same physical mapping, serial and 2-rank, with restart on and off; 2-rank runs also swap rank ownership; a multi-target fixture (`-t05-multi`) covers one target face overlapping two source faces; `T05-unmatched-source-*` checks that a required unmatched source face fails on all ranks, with restart on and off |
| T06 | Fully covered `2.5e-7` fine face in both ownership directions; no target face left uncovered |
| T07 | Normal-gap discovery (gap 0.03, `maxDistance2 0.01`): coverage and pair-set equality on 2 and 4 ranks; split-grid pair-set oracle on 2 ranks (4 ranks deferred) |
| T08-T09 | Explicitly not run in S02; full topology, state, FV, and GAMG matrices are scheduled for S03 |
| T10 | Partial-coverage diagnostics |
| T11 | Scale/offset matrix scheduled for S03 and hardening in S10 |
| T12 | Empty-target input |
| T13 | Piecewise-constant area-average control |
| T14-T19 | `pureZoneMixture` dictionaries, material properties, cell-zone coverage, analytical references, and static meshes validated; `solidFoam` is not run until S03 |

The T05 case checks seed selection on the two-island fixture (`-t05-order`).
Each source island must be seeded from a target face it overlaps, and the four
face orderings must give the same physical mapping (`T05_PHYSICAL_MAPPING`).
That check uses stable face IDs and runs serially and under MPI: each rank
reports its local source-face weights, and the result is reduced over ranks.
The records run with restart off (`-no-restart`, the S02 baseline), with restart
on (serial, 2-rank distributed, and 2-rank swapped ownership), and on a
multi-target fixture (`-t05-multi`). Restart on the two-island fixture does not
exercise restart recovery, because every face is already fully covered; the
restart-recovery cases are T01-T04, T06-T08, T10, T11, and T13.

The negative case (`-unmatched-source`) uses `requireMatch true` with one source
face that has no overlap. It must fail on every rank, with restart on and off.
With restart on, an unmapped required face is a collective fatal after the
restart pass. This is the S07 behaviour (decision D1). Previously, with restart
on, the face stayed unmapped without error.

The T07 case checks normal-gap discovery (gap `0.03`, `-maxDistance2 0.01`,
`requireMatch true`). A pair accepted by the serial walk must also be accepted
under MPI. The broad phase in `calcProcMap()` enlarges the target-face box by a
distance bound derived from `maxDistance2` and `minCosAngle` when both are set.
Otherwise it sends every non-empty target face to every source processor. The T07
assertions are coverage (`T07_CALCULATED`), pair-set equality against the
all-pairs oracle (`T07_PAIR_MISMATCH`), and pair count (`T07_PAIR_COUNT`).
Options: `-tilt <deg>` tilts the target plane, `-minCosAngle <cos>` sets the
normal-cosine limit, and `-split-grid N` splits an `N x N` source patch across
ranks (T07 only).

Gated cases: fallback and bounded modes on 2 and 4 ranks (coverage cases), and
split-grid fallback and bounded modes on 2 ranks. Split-grid on 4 ranks is not
gated. It shows 1560 per-face mismatches on receive-only ranks. The original
AABB `calcProcMap()` shows the same defect, so it is pre-existing and not
caused by the bound change. It is deferred with a follow-up. Traffic (send-map
entries, 100x100 split grid, gap 0.03): fallback 20000 (2 ranks) and 40000
(4 ranks); bounded 12700 and 18100.

The bound is derived from the kernel test, not proven for oblique source
normals. The fallback covers any case the bound does not. It does not replace the broader all-pairs discovery oracle and rank-topology
matrix in the approved plan. This is the targeted S02 baseline gate, not the
full test matrix in every S02 substep. T02 currently checks the midpoint
fraction; T05 checks seed selection with and without restart; T06 checks the fully covered fine face; T04
checks vector fallback values; and T12 checks an empty target. Their remaining
fraction, ordering, cutoff, field-type, and invalid-input variants are not
counted as passed. Likewise, T13 checks piecewise-constant area averaging. It
does not claim second-order accuracy.

T20 is the S08 symmetric overlap-retention negative control. A unit source face
and a target face of area 0.75 overlap in a strip of width `w` and length 0.5,
so the overlap is `0.5w`. The retention threshold is `1e-6 * min(area)` = `7.5e-7`.
With `-sliver-width 1.25e-6` (overlap `6.25e-7`) the overlap must be rejected in
both ownership directions (`-swap`), so both coverages are zero. With
`-sliver-width 2e-6` (overlap `1e-6`) it must be retained, and the unit and
sliver coverages are checked against their expected values. T20 is run only
with `--tests T20`; it is not part of the S02 gate. `-grid N` (T01 only) builds
a uniform N x N grid for stencil-size and mapping-timing comparison; it
returns after printing the stencil counts and does not run the generic T01
coverage checks. `-grid-small s` (with `-grid N`, T01 only) splits each cell
3x3 with centre faces of side `s` relative to the cell, giving target faces
far smaller than their source faces. `-grid-sliver w` (with `-grid N`, T01
only) adds thin target strips of relative width `w`. Both print an
`AMI_STENCIL` line with link counts and uncovered source and target faces.
Slivers below about `1e-6` of the cell size are removed by the clipping snap
in `faceAreaIntersect`, so `-grid-sliver` cannot test narrower slivers.

## Solid fixture contract

`cases/solidFoam/variants.json` records each variant's material properties,
mesh resolution, boundary temperatures, and reference dimensions. The runner
checks the rendered `pureZoneMixture` dictionary with `foamDictionary`, checks
that `left` and `right` cell zones cover each cell exactly once, runs
`checkMesh`, and records the analytical reference metrics from `solidMetrics.py`.
The fixture template has no unzoned cells. Its `none` material entry is present
for dictionary completeness and is not used by any cell.

`test_solidMetrics.py` checks the independent references. T15 uses the
series-resistance heat rate and junction temperature. T16 records the
volumetric heat capacities. T17 records the full rotated conductivity tensor
and flux with an off-diagonal contribution. The accepted solid error limits
are recorded with those references; transient, anisotropic, and moving
accuracy limits remain subject to the C0 characterization required by the
plan.

S02 validates only inputs and static meshes for T14-T19. It does not prove
runtime selection of `pureZoneMixture`, run a thermal solve, or validate moving
mesh/topology updates. Those checks remain explicit S03 work. Do not use the
fixture-only status as a solver result.

## Files

- `Test-AMIInterpolation.C` contains the named OpenFOAM assertions.
- `runTests.py` builds and runs the driver, checks C0 provenance, handles
  target provenance, handles timeouts, validates solid fixtures, and writes
  machine-readable results for S02 and S03.
- `runS03ApplicationScreens.py` runs isolated B0/C0 solid and copied-tutorial
  screens and writes application-level metrics and input hashes.
- `solidMetrics.py` contains independent solid reference calculations and
  acceptance limits.
- `test_solidMetrics.py` tests those calculations with Python's standard
  `unittest` module.
- `cases/` contains canonical dictionaries and the parameterized solid case.
  Generated meshes, processor data, executables, and logs belong in the
  external AMI project workspace, not in this directory.
