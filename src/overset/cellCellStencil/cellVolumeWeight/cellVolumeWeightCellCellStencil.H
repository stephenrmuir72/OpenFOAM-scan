/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::cellCellStencils::cellVolumeWeight

Description
    Volume-overlap interpolation stencil. Default positive weights preserve
    constants and boundedness, but are not generally affine-exact or a globally
    conservative overset coupling.

    Optional controls:
      donorSelection       legacy (default) | coverage
      donorFallback        legacy (default) | alternateZones | sameZone
      minCoverage          [0,1], default 0
      minRetainedFraction  [0,1], default 0
      weightingScheme      volume (default) | linearExact
      allowNegativeWeights false (default)
      stencilExpansionRings 0 | 1 (default)
      layerDonorMetric      legacy (default) | overlapWeighted
      processorMapMethod   AABB (default) | LOD

    Robustness controls (safe defaults):
      checkAnchoring       true (default): detect closed interpolation cycles
                           (fully interpolated cells without a positive-weight
                           donor path to a solved cell; singular I - W block)
                           and repair them by dropping unanchored donors,
                           re-selecting a candidate or demoting the cell.
                           Also reverts affine stencils that close a cycle.
      checkPockets         none | report (default) | fix: solved regions with
                           neither an interpolated cell nor a physical patch
                           (pure Neumann, singular pressure). fix blanks them.
      patchCutting         legacy (default) | overlapFraction | coverage.
                           legacy blanks any acceptor overlapping a donor cell
                           next to a physical patch (removes all acceptors in
                           narrow gaps, e.g. press/table). overlapFraction
                           blanks if the PATCH-donor overlap exceeds
                           patchOverlapFraction*V; coverage blanks only cells
                           not covered by the donor zone (1-overlapTolerance).
      patchOverlapFraction [0,1], default 0
      patchCutHysteresis   (0,1], default 1 (none): threshold factor for
                           cells that were holes, against hole flicker.
      holePoints           ((zone (x y z)) ...): points inside solid bodies in
                           the points0 (reference) configuration of the zone.
                           The rigidly fitted (Kabsch) current position seeds
                           forced holes in the other zones, also for bodies
                           crossing the domain boundary.
      stopIneligibleHoleGrowth false (default): front cells without donors
                           become porous (solved) instead of growing the hole.
      donorPriority        (z0 z1 ...): donor zones in decreasing priority;
                           default () = nearest zone index.
      donorSwitchMargin    >= 0, default 0: coverage margin favouring the
                           previous donor zone (hysteresis).

    linearExact controls:
      affineCorrection     minNorm (default) | weighted | positive. weighted
                           minimises sum dw^2/w0 (keeps small donors small);
                           positive also clips negative weights and
                           re-projects (bounded, may fall back).
      maxWeightAmplification  > 0, default 1/sqrt(eps): bound on sum |w|
      affineSvdTolerance   (0,1), default 100 eps
      affineScaling        perAxis (default) | isotropic
    allowNegativeWeights must stay false for bounded fields (e.g. VOF alpha).

    Performance controls (exact unless noted):
      cacheStaticPairs     false (default): reuse overlap data of zone pairs
                           whose points did not move (bitwise), e.g. table and
                           background while only the press moves.
      skipStaticUpdate     false (default): keep the previous stencil when no
                           zone moved.
      pairOrientation      zoneOrder (default) | smallerTarget (search tree on
                           the smaller zone; equal up to roundoff).
      overlapShortcut      false (default): overlap of a cell strictly inside a
                           convex cell is its volume (equal up to roundoff).

    Accuracy: volume weights are first-order (piecewise constant). Donor
    cells much coarser than acceptors smear gradients; keep donor/acceptor
    size ratios near 1 in the overlap and use linearExact for second-order
    interpolation. For moving bodies near walls use at least two donor cell
    layers between overset boundaries and walls, overlapTolerance around
    1e-6..1e-3 and, for coverage selection, minCoverage 0.9..0.99.

    The affine option falls back to validated overlap weights when its moment
    constraints or boundedness policy cannot be satisfied. Recovery uses IDW
    as its fallback instead of geometric volumes. Signed weights require
    explicit opt-in and can produce overshoots or negative interpolated fields.

    Coverage uses the overlap engine's first-vertex fan volume, which can differ
    from FV volume on warped faces. FV volumes remain unchanged. Inverted fan
    tets (concave or warped faces) are handled by a signed decomposition in
    both the volume and the overlap engine; flat fan tets are ignored.
    Optional coverage gates apply before fringe selection and after filtering.
    IDW has separate provenance and is never interpreted as geometric coverage.

    Finite restart classifications are rounded/clamped with global warnings;
    non-finite and incomplete non-empty-rank restart data are errors.

SourceFiles
    cellVolumeWeightCellCellStencil.C

\*---------------------------------------------------------------------------*/

#ifndef cellCellStencils_cellVolumeWeight_H
#define cellCellStencils_cellVolumeWeight_H

#include "cellCellStencil.H"
#include "volFields.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMeshSubset;

namespace cellCellStencils
{

/*---------------------------------------------------------------------------*\
                     Class cellVolumeWeight Declaration
\*---------------------------------------------------------------------------*/

class cellVolumeWeight
:
    public cellCellStencil
{
protected:

    // Static data members

        //- Default overlap tolerance. Fraction of volume
        static scalar defaultOverlapTolerance_;


    // Protected data

        //- Dictionary of motion control parameters
        const dictionary dict_;

        //- Tolerance for volume overlap. Fraction of volume
        scalar overlapTolerance_;

        //- Per cell the cell type
        labelList cellTypes_;

        //- Indices of interpolated cells
        labelList interpolationCells_;

        //- Fetch interpolated cells
        autoPtr<mapDistribute> cellInterpolationMap_;

        //- Interpolation stencil
        labelListList cellStencil_;

        //- Interpolation weights
        List<scalarList> cellInterpolationWeights_;

        //- Amount of interpolation
        volScalarField cellInterpolationWeight_;

        //- Allow interpolated as donors
        const bool allowInterpolatedDonors_;

        const word donorSelection_;
        const word donorFallback_;
        const scalar minCoverage_;
        const scalar minRetainedFraction_;
        const word weightingScheme_;
        const bool allowNegativeWeights_;
        const label stencilExpansionRings_;
        const word layerDonorMetric_;
        const word processorMapMethod_;

        //- Detect and repair closed interpolation cycles (default true)
        const bool checkAnchoring_;

        //- Isolated solved-pocket check: none | report (default) | fix
        const word checkPockets_;

        //- Acceptor cut next to donor walls:
        //  legacy (default) | overlapFraction | coverage
        const word patchCutting_;

        //- overlapFraction: cut if PATCH-donor overlap > fraction*V
        const scalar patchOverlapFraction_;

        //- Threshold factor in (0,1] for cells that were holes (1 = none)
        const scalar patchCutHysteresis_;

        //- Optional points inside solid bodies, (zone point) in the
        //  reference (points0) configuration of the zone. Cells of other
        //  zones containing the rigidly moved point seed forced holes.
        List<Tuple2<label, point>> holePoints_;

        //- Reference points for holePoints (lazily read)
        autoPtr<pointField> holePoints0_;

        //- Front cells without donors become POROUS instead of growing the
        //  hole (default false = legacy)
        const bool stopIneligibleHoleGrowth_;

        //- Donor zones in decreasing priority (default empty = legacy
        //  nearest-zone-index rule)
        const labelList donorPriority_;

        //- Coverage margin favouring the previous donor zone (default 0)
        const scalar donorSwitchMargin_;

        //- Donor zone of interpolated cells from the previous update
        labelList donorZone_;

        //- linearExact correction: minNorm (default) | weighted | positive
        const word affineCorrection_;

        //- Maximum sum |w| of corrected weights (default 1/sqrt(eps))
        const scalar maxWeightAmplification_;

        //- Relative SVD truncation (default 100 eps)
        const scalar affineSvdTolerance_;

        //- Coordinate scaling: perAxis (default) | isotropic
        const word affineScaling_;

        //- Overlap-tessellation cell volumes and roundoff (lazily evaluated
        //  after the first update; -1 = not yet evaluated)
        mutable scalarField geometricVolumes_;
        mutable scalarField geometricRoundoff_;
        mutable labelList weightFallbackCounts_;
        mutable bool weightFallbackReported_;
        mutable bool ambiguityReported_;
        bool geometryMismatchReported_;
        bool thinMeshReported_;

        struct donorCandidate
        {
            label zone;
            labelList donors;
            scalarList volumes;
        };

        //- Overlap data of one side of a zone pair (acceptor side)
        struct pairSide
        {
            labelListList addressing;
            scalarListList weights;
            labelList donors;
            labelList patchTypes;
            scalarList patchOverlap;
        };

        //- Reuse overlap data of zone pairs whose zones did not move
        //  (default false)
        const bool cacheStaticPairs_;

        //- Return without recomputation if no zone moved (default false)
        const bool skipStaticUpdate_;

        //- zoneOrder (default) | smallerTarget
        const word pairOrientation_;

        //- Containment shortcut in the overlap computation (default false)
        const bool overlapShortcut_;

        //- Cached zone-pair data (index srcZone*nZones + tgtZone)
        PtrList<FixedList<pairSide, 2>> pairCache_;

        //- Zone of each point (-1 none, -2 several) for motion detection
        labelList pointZone_;

        //- Points at the previous update for motion detection
        pointField lastPoints_;

        //- A stencil has been computed by update()
        bool haveStencil_;

   // Protected Member Functions

        static bool normaliseVolumeWeights(scalarList& weights);

        static label decodeRestartType
        (
            const scalar value,
            bool& rounded,
            bool& clamped
        );

        void refreshGeometry();

        //- Evaluate the geometric volume of a cell
        void evaluateGeometry(const label celli) const;

        //- Geometric (overlap tessellation) volume of a cell
        scalar geometricVolume(const label celli) const
        {
            if (geometricVolumes_[celli] < 0) evaluateGeometry(celli);
            return geometricVolumes_[celli];
        }

        //- Relative roundoff allowance of the geometric volume
        scalar geometricRoundoff(const label celli) const
        {
            if (geometricRoundoff_[celli] < 0) evaluateGeometry(celli);
            return geometricRoundoff_[celli];
        }

        bool candidateCoverage
        (
            const label celli,
            const scalarList& original,
            const scalarList& retained,
            scalar& coverage
        ) const;

        bool betterCandidate
        (
            const label celli,
            const label zone,
            const scalar coverage,
            const label bestZone,
            const scalar bestCoverage
        ) const;

        //- Zone preference (donorPriority or legacy nearest zone index)
        bool preferZone
        (
            const label acceptorZone,
            const label zone,
            const label currentZone
        ) const;

        //- Previous donor zone of a cell (-1 if unknown)
        label previousDonorZone(const label celli) const
        {
            return donorZone_.size() == mesh_.nCells() ? donorZone_[celli] : -1;
        }

        scalar weightedDonorVolume
        (
            const pointList& centres,
            const scalarList& volumes,
            scalarList weights
        ) const;

        static void distanceWeights
        (
            const point& sample,
            const pointList& centres,
            scalarList& weights
        );

        bool affineWeights
        (
            const point& sample,
            const pointList& centres,
            scalarList& weights
        ) const;

        void reportWeightFallbacks
        (
            label nAffineExact,
            label nAffineFallback
        ) const;

        void finaliseDonors
        (
            const globalIndex& globalCells,
            const boolList& previouslyBlanked,
            const boolList& sameZoneFallback,
            const List<DynamicList<donorCandidate>>& candidates,
            labelList& types,
            labelListList& stencil,
            scalarListList& weights,
            labelList& donorZones,
            const boolList* donorExclusions = nullptr
        ) const;

        //- Per cell whether its value is (transitively) determined by a
        //  solved cell: non-interpolated, partially blended, or
        //  interpolated from at least one anchored donor
        boolList anchoredCells
        (
            const globalIndex& globalCells,
            const labelList& types,
            const labelListList& stencil,
            const scalarListList& weights,
            const scalarField& blend
        ) const;

        //- Repair unanchored interpolated cells (closed interpolation
        //  cycles): drop unanchored donors, else re-select an anchored
        //  candidate, else demote. Returns true if any cell was demoted.
        //  counts: dropped-donor repairs, re-selections, demotions
        bool repairUnanchored
        (
            const globalIndex& globalCells,
            const boolList& anchored,
            const boolList& previouslyBlanked,
            const boolList& sameZoneFallback,
            const List<DynamicList<donorCandidate>>& candidates,
            const bool forceDemote,
            labelList& types,
            labelListList& stencil,
            scalarListList& weights,
            scalarField& blend,
            labelList& donorZones,
            labelList& counts
        ) const;

        //- Mark solved cells in connected regions (faces blocked at
        //  HOLE|non-HOLE) that have neither an interpolated cell with
        //  positive blend nor a real (non-overset, non-constraint) patch
        //  face. Returns the global number of such cells.
        label isolatedPockets
        (
            const labelList& types,
            const scalarField& blend,
            boolList& isolated
        ) const;

        //- Find cells next to cells of type PATCH
        void findHoles
        (
            const globalIndex& globalCells,
            const fvMesh& mesh,
            const labelList& zoneID,
            const labelListList& stencil,
            labelList& cellTypes,
            const boolList& holeSeeds = boolList()
        ) const;

        //- Local copy of cellCellStencil::walkFront (identical results with
        //  stopIneligibleHoleGrowth false). Front iteration visits set bits
        //  only. With stopIneligibleHoleGrowth, cells without donors reached
        //  by the front become POROUS and do not re-seed the front.
        //  Returns the number of such cells (holes or porous) in the
        //  selected layer.
        label walkFrontLocal
        (
            const globalIndex& globalCells,
            const scalar layerRelax,
            const labelListList& allStencil,
            labelList& allCellTypes,
            scalarField& allWeight,
            const scalarList& compactCellVol,
            const labelListList& compactStencil,
            const labelList& zoneID,
            const label holeLayers,
            const label useLayer,
            boolList* protectedHoles = nullptr,
            const scalarList* donorVolumeRatio = nullptr
        ) const;

        //- Cells of other zones containing the (rigidly moved) holePoints
        void holePointCells
        (
            const PtrList<fvMeshSubset>& meshParts,
            const labelList& zoneID,
            boolList& seeds
        );

        //- according to additionalDocumentation/MEJ_oversetMesh.txt
        void markPatchCells
        (
            const fvMesh& mesh,
            const labelList& cellMap,
            labelList& patchCellTypes
        ) const;

        void combineCellTypes
        (
            const label subZoneID,
            const fvMesh& subMesh,
            const labelList& subCellMap,

            const label donorZoneID,
            const labelListList& toOtherCells,
            const List<scalarList>& weights,
            const labelList& otherCells,
            const labelList& interpolatedOtherPatchTypes,

            labelListList& allStencil,
            scalarListList& allWeights,
            labelList& allCellTypes,
            labelList& allDonorID,
            List<DynamicList<donorCandidate>>* candidates = nullptr,
            const scalarList& patchOverlap = scalarList()
        ) const;

        //- Per-acceptor overlap volume with donor cells of type PATCH
        //  (used by patchCutting overlapFraction)
        void patchOverlapVolumes
        (
            const autoPtr<mapDistribute>& mapPtr,
            const labelListList& addressing,
            const scalarListList& weights,
            const labelList& patchTypes,
            scalarList& result
        ) const;

        //- interpolate (= combine) patch types
        void interpolatePatchTypes
        (
            const labelListList& addressing,
            const labelList& patchTypes,
            labelList& result
        ) const;

        //- interpolate (= combine) patch types
        void interpolatePatchTypes
        (
            const autoPtr<mapDistribute>& mapPtr,
            const labelListList& addressing,
            const labelList& patchTypes,
            labelList& result
        ) const;


private:

    // Private Member Functions

        //- No copy construct
        cellVolumeWeight(const cellVolumeWeight&) = delete;

        //- No copy assignment
        void operator=(const cellVolumeWeight&) = delete;


public:

    //- Runtime type information
    TypeName("cellVolumeWeight");


    // Constructors

        //- Construct from fvMesh
        cellVolumeWeight(const fvMesh&, const dictionary&, const bool doUpdate);


    //- Destructor
    virtual ~cellVolumeWeight();


    // Member Functions

        //- Access to volume overlap tolerance
        scalar overlapTolerance() const
        {
            return overlapTolerance_;
        }

        //- Update stencils. Return false if nothing changed.
        virtual bool update();

        //- Return the cell type list
        virtual const labelUList& cellTypes() const
        {
            return cellTypes_;
        }

        //- Indices of interpolated cells
        virtual const labelUList& interpolationCells() const
        {
            return interpolationCells_;
        }

        //- Return a communication schedule
        virtual const mapDistribute& cellInterpolationMap() const
        {
            if (!cellInterpolationMap_)
            {
                const_cast<cellVolumeWeight&>(*this).update();
            }
            return *cellInterpolationMap_;
        }

        //- Per interpolated cell the neighbour cells (in terms of slots as
        //  constructed by above cellInterpolationMap) to interpolate
        virtual const labelListList& cellStencil() const
        {
            return cellStencil_;
        }

        //- Weights for cellStencil
        virtual const List<scalarList>& cellInterpolationWeights() const
        {
            return cellInterpolationWeights_;
        }

        //- Per interpolated cell the interpolation factor. (0 = use
        //  calculated, 1 = use interpolated)
        virtual const scalarList& cellInterpolationWeight() const
        {
            return cellInterpolationWeight_;
        }

        //- Same-zone recovery weights (IDW baseline, optional affine correction)
        virtual void stencilWeights
        (
            const point& sample,
            const pointList& donorCcs,
            scalarList& weights
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace cellCellStencils
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
