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

    The affine option falls back to validated overlap weights when its moment
    constraints or boundedness policy cannot be satisfied. Recovery uses IDW
    as its fallback instead of geometric volumes. Signed weights require
    explicit opt-in and can produce overshoots or negative interpolated fields.

    Coverage uses the overlap engine's first-vertex fan volume, which can differ
    from FV volume on warped faces. FV volumes remain unchanged. Positive
    volume-carrying fan tets are required; exactly collinear face triangles are
    ignored. Numerically unresolved tetrahedron areas are errors.
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
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
        scalarField geometricVolumes_;
        scalarField geometricRoundoff_;
        mutable labelList weightFallbackCounts_;
        mutable bool weightFallbackReported_;

        struct donorCandidate
        {
            label zone;
            labelList donors;
            scalarList volumes;
        };

   // Protected Member Functions

        static bool normaliseVolumeWeights(scalarList& weights);

        static label decodeRestartType
        (
            const scalar value,
            bool& rounded,
            bool& clamped
        );

        void refreshGeometry();

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

        static scalar weightedDonorVolume
        (
            const pointList& centres,
            const scalarList& volumes,
            scalarList weights
        );

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

        //- Find cells next to cells of type PATCH
        void findHoles
        (
            const globalIndex& globalCells,
            const fvMesh& mesh,
            const labelList& zoneID,
            const labelListList& stencil,
            labelList& cellTypes
        ) const;

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
            List<DynamicList<donorCandidate>>* candidates = nullptr
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
