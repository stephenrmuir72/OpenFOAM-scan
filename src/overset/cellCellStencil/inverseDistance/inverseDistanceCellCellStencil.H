/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2022 OpenCFD Ltd.
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
    Foam::cellCellStencils::inverseDistance

Description
    Inverse-distance-weighted interpolation stencil.

    The optional weightingScheme selects the weight calculation. The default
    is inverse-distance weighting; least-squares weighting can reproduce
    affine fields when the donor geometry is sufficiently well-conditioned.
    Degenerate stencils use nearest-donor injection and are not generally
    affine-exact.

    hole finding:
    - mark boundary faces on helper (voxel) mesh
    - mark any cell overlaying these voxels
    - use flood filling to find any unreachable cell
    Alternative is to use an octree of the boundary faces and determine
    directly for all cells whether we are outside. Might be slow though.

SourceFiles
    inverseDistanceCellCellStencil.C

\*---------------------------------------------------------------------------*/

#ifndef cellCellStencils_inverseDistance_H
#define cellCellStencils_inverseDistance_H

#include "cellCellStencil.H"
#include "volFields.H"
#include "labelVector.H"
#include "treeBoundBox.H"
#include "pointList.H"
#include "globalIndex.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMeshSubset;

namespace patchFaceIntersection
{
class triangleSearch;
}

namespace cellCellStencils
{

/*---------------------------------------------------------------------------*\
                     Class inverseDistance Declaration
\*---------------------------------------------------------------------------*/

class inverseDistance
:
    public cellCellStencil
{
    // Private Member Functions

        //- No copy construct
        inverseDistance(const inverseDistance&) = delete;

        //- No copy assignment
        void operator=(const inverseDistance&) = delete;


protected:

    // Protected data

        //- Dictionary of motion control parameters
        const dictionary dict_;

        enum class weightingScheme
        {
            INVERSE_DISTANCE,
            LINEAR_EXACT
        };
        static const Enum<weightingScheme> weightingSchemeNames_;
        const weightingScheme weightingScheme_;
        mutable labelList weightFallbackCounts_;
        mutable bool weightFallbackReported_;

        //- Allow holes as donors
        const bool allowHoleDonors_;

        //- Allow interpolared as donors
        const bool allowInterpolatedDonors_;

        //- Geometry tolerance scaled by the current mesh bounds span.
        //  The 'tolerance' dictionary entry defaults to 1e-10; zero requests
        //  exact predicates. Refreshed at each geometry update.
        vector smallVec_;

        //- Per cell the cell type
        labelList cellTypes_;

        //- Indices of interpolated cells
        labelList interpolationCells_;

        //- Fetch interpolated cells
        autoPtr<mapDistribute> cellInterpolationMap_;

        //- Interpolation stencil
        labelListList cellStencil_;

        //- Interpolation weights
        scalarListList cellInterpolationWeights_;

        //- Amount of interpolation
        volScalarField cellInterpolationWeight_;

        //- Data for a set of interpolated/donor set
        struct interpolatedDonorSet
        {
            scalar averVolRatio;
            labelList cellTypes;
        };

        //- Candidate donor cell and the mesh zone it belongs to
        struct donorCandidate
        {
            label zone;
            label globalCell;

            donorCandidate()
            :
                zone(-1),
                globalCell(-1)
            {}

            donorCandidate(const label zonei, const label globalCelli)
            :
                zone(zonei),
                globalCell(globalCelli)
            {}
        };


   // Protected Member Functions

        void nearestWeights
        (
            const point& sample,
            const pointList& donorCcs,
            scalarList& weights,
            const label reason
        ) const;

        void reportWeightFallbacks() const;

        void restoreReactivationCandidates
        (
            const globalIndex& globalCells,
            const labelList& zoneID,
            const labelList& patchTypes,
            const boolList& previouslyBlanked,
            const boolList& candidates,
            boolList& sameZoneFallback,
            const labelListList& reactivationStencil,
            const boolList& protectedHoles,
            labelList& cellTypes,
            labelListList& stencil,
            scalarField& interpolationWeight
        ) const;

        //- Refresh the geometry tolerance derived from the current mesh span
        void updateGeometryTolerance();

        //- Check whether a source cell is a valid donor for an acceptor
        bool donorEligible
        (
            const label donorType,
            const label donorPatchType,
            const label acceptorPatchType
        ) const;

        //- Retain a donor candidate and update the best donor zone
        void recordDonorCandidate
        (
            const label acceptorCell,
            const label acceptorZone,
            const label donorZone,
            const label globalDonor,
            List<DynamicList<donorCandidate>>& allCandidates,
            labelList& allDonor
        ) const;

        //- Rank donor candidates using the existing zone preference
        void orderDonorCandidates
        (
            const labelList& zoneID,
            List<DynamicList<donorCandidate>>& allCandidates,
            labelListList& allStencil,
            labelList& allDonor
        ) const;

        //- Re-run hole classification and front walking when a donor fails
        void classifyAndWalkFront
        (
            const globalIndex& globalCells,
            const labelList& zoneID,
            const labelList& allPatchTypes,
            const boolList& previouslyBlanked,
            const List<DynamicList<donorCandidate>>& allCandidates,
            const scalar layerRelax,
            const label holeLayers,
            const label useLayer,
            labelList& allCellTypes,
            labelListList& allStencil,
            scalarField& allWeight,
            boolList& reactivationCandidates,
            boolList& sameZoneFallback,
            labelListList& reactivationStencil,
            boolList& protectedHoles
        ) const;

        // Voxel representation

            //- Convert ijk indices into single index
            static label index(const labelVector& nDivs, const labelVector&);

            //- Convert single index into ijk
            static labelVector index3(const labelVector& nDivs, const label);

            //- Convert coordinate into ijk
            static labelVector index3
            (
                const boundBox& bb,
                const labelVector& nDivs,
                const point& pt
            );

            //- Convert index back into coordinate
            static point position
            (
                const boundBox& bb,
                const labelVector& nDivs,
                const label boxI
            );

            //- Fill all elements overlapping subBb with value val
            static void fill
            (
                PackedList<2>& elems,
                const boundBox& bb,
                const labelVector& nDivs,
                const boundBox& subBb,
                const unsigned int val
            );

            //- Is any voxel inside subBb set to val
            static bool overlaps
            (
                const boundBox& bb,
                const labelVector& nDivs,
                const PackedList<2>& voxels,
                const treeBoundBox& subBb,
                const unsigned int val
            );

        //- Mark voxels of patchTypes with type of patch face
        static void markBoundaries
        (
            const fvMesh& mesh,
            const vector& smallVec,

            const boundBox& bb,
            const labelVector& nDivs,
            PackedList<2>& patchTypes,
            const labelList& cellMap,
            labelList& patchCellTypes
        );

        //- Mark all cells overlapping (a voxel covered by) a src patch
        //  with type HOLE
        void markPatchesAsHoles
        (
            PstreamBuffers& pBufs,

            // Mesh bb and data
            const PtrList<fvMeshSubset>& meshParts,

            // Helper mesh for patches
            const List<treeBoundBoxList>& patchBb,
            const List<labelVector>& patchDivisions,
            const PtrList<PackedList<2>>& patchParts,
            const patchFaceIntersection::triangleSearch& srcSearch,

            const label srcI,
            const label tgtI,
            labelList& allCellTypes
        ) const;

        //- If multiple donors meshes: decide which is best
        bool betterDonor
        (
            const label destMesh,
            const label currentDonorMesh,
            const label newDonorMesh
        ) const;

        //- Determine donors for all tgt cells
        void markDonors
        (
            const globalIndex& globalCells,
            PstreamBuffers& pBufs,
            const PtrList<fvMeshSubset>& meshParts,
            const List<treeBoundBoxList>& meshBb,
            const labelList& allCellTypes,
            const labelList& allPatchTypes,

            const label srcI,
            const label tgtI,
            List<DynamicList<donorCandidate>>& allCandidates,
            labelList& allDonor
        ) const;

        //- Replacement of regionSplit
        void uncompactedRegionSplit
        (
            const fvMesh& mesh,
            const globalIndex& globalFaces,
            const label nZones,
            const labelList& zoneID,
            const labelList& cellTypes,
            const boolList& isBlockedFace,
            labelList& cellRegion
        ) const;

        autoPtr<globalIndex> compactedRegionSplit
        (
            const fvMesh& mesh,
            const globalIndex& globalFaces,
            labelList& cellRegion
        ) const;

        //- Do flood filling to detect unreachable (from patches) sections
        //  of mesh
        void findHoles
        (
            const globalIndex& globalCells,
            const fvMesh& mesh,
            const labelList& zoneID,
            const labelListList& stencil,
            labelList& cellTypes
        ) const;

        //- Seed faces of cell with wantedFraction (if higher than current)
        void seedCell
        (
            const label cellI,
            const scalar wantedFraction,
            bitSet& isFront,
            scalarField& fraction
        ) const;

        //- Create stencil starting from the donor containing the acceptor
        //- allowHoleDonors : allow donors of type HOLE (otherwise are filtered
        //- out)
        virtual void createStencil
        (
            const globalIndex&,
            const bool allowHoleDonors
        );

        //- Make holes next to live ones type SPECIAL and include in
        //- interpolation stencil
        void holeExtrapolationStencil(const globalIndex& globalCells);

        //- Calculate least-squares weights for a single acceptor
        void linearExactWeights
        (
            const point& sample,
            const pointList& donorCcs,
            scalarList& weights
        ) const;


public:

    //- Runtime type information
    TypeName("inverseDistance");


    // Constructors

        //- Construct from fvMesh
        inverseDistance(const fvMesh&, const dictionary&, const bool);


    //- Destructor
    virtual ~inverseDistance();


    // Member Functions

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
                const_cast<inverseDistance&>(*this).update();
            }
            return cellInterpolationMap_();
        }

        //- Per interpolated cell the neighbour cells (in terms of slots as
        //  constructed by above cellInterpolationMap) to interpolate
        virtual const labelListList& cellStencil() const
        {
            return cellStencil_;
        }

        //- Weights for cellStencil
        virtual const scalarListList& cellInterpolationWeights() const
        {
            return cellInterpolationWeights_;
        }

        //- Per interpolated cell the interpolation factor. (0 = use
        //  calculated, 1 = use interpolated)
        virtual const scalarList& cellInterpolationWeight() const
        {
            return cellInterpolationWeight_;
        }

        //- Calculate inverse distance weights for a single acceptor
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
