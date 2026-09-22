/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Patrick Wilhelm, Data Lab Hell
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
    Foam::multiRegionCoupledDecomp

Description
    Single globally-consistent decomposition for implicitly coupled
    multi-region cases, selectable as \c multiRegionCoupled.

    Welds the cells either side of every coupled (mapped) face onto one
    processor, so the assembled multi-region matrix
    (lduPrimitiveMeshAssembly) keeps each such face local, then balances
    the welded super-cells with any sub-method. \c implicitPatches and
    \c explicitPatches filter which mapped patches are welded;
    \c regionWeights biases the balance towards solver work.

    Serial only. \c decomposePar \c -allRegions places every region of
    \c regionProperties in one pass, from a single decomposition built from
    the first region decomposed; a later region with different coefficients,
    or one selecting another method, is refused. Requires conformal
    \c nearestPatchFace coupling.
    \c decomposeParDict \c constraints are not supported.

    \c patchWeights adds a per-face cost, nested per region since patch
    names are unique only within a region.

    Method coefficients:
    \table
        Property        | Description                       | Required | Default
        method          | sub-method balancing the welded graph | yes | -
        implicitPatches | weld only patches matching these names <!--
                        --> | no | all mapped
        explicitPatches | never weld patches matching these names | no | none
        regionWeights   | per-region cost factors (sub-dictionary) <!--
                        --> | no | all 1
        patchWeights    | per-face cost factors, a sub-dictionary of <!--
                        --> patch weights per region | no | none
    \endtable

SourceFiles
    multiRegionCoupledDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_multiRegionCoupledDecomp_H
#define Foam_multiRegionCoupledDecomp_H

#include "decompositionMethod.H"
#include "regIOobject.H"
#include "UPtrList.H"
#include "wordRes.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class multiRegionCoupledDecomp Declaration
\*---------------------------------------------------------------------------*/

class multiRegionCoupledDecomp
:
    public decompositionMethod
{
    // Private Classes

        //- Welded decomposition over all regions, shared by
        //  per-region decompose() calls so registered on Time
        class weldedDecomposition
        :
            public regIOobject
        {
        public:

            //- Region names, in regionProperties order
            wordList regions;

            //- Start of each region in the global cell numbering
            labelList regionOffsets;

            //- Processor per global cell
            labelList decomp;

            //- Region whose decompose() built this
            word ownerName;

            //- Coefficients it was built with; later regions must match
            dictionary coeffs;

            TypeNameNoDebug("multiRegionCoupledWeld");

            explicit weldedDecomposition(const objectRegistry& obr);

            //- Write data is no-op
            virtual bool writeData(Ostream&) const
            {
                return true;
            }
        };


    // Private Data

        //- Coefficients dictionary for the sub-method
        dictionary methodDict_;

        //- Sub-method that balances the welded super-cells
        autoPtr<decompositionMethod> method_;

        //- Optional weld filters, tested on either side of an interface
        wordRes implicitPatches_;
        wordRes explicitPatches_;

        //- Optional per-region cost factors (region name -> weight)
        dictionary regionWeights_;

        //- Optional per-face cost factors
        dictionary patchWeights_;


    // Private Member Functions

        //- The welded decomposition, built on first use
        const weldedDecomposition& weldedDecomp(const polyMesh& mesh) const;

        //- Build the welded decomposition over all regions, reading the
        //  region meshes
        void buildGlobalDecomposition
        (
            const polyMesh& mesh,
            weldedDecomposition& weld
        ) const;

        //- Refuse a region whose own decomposeParDict selects another
        //  method
        void checkRegionMethods
        (
            const polyMesh& mesh,
            const UPtrList<const polyMesh>& meshes
        ) const;


public:

    // Generated Methods

        //- No copy construct
        multiRegionCoupledDecomp(const multiRegionCoupledDecomp&) = delete;

        //- No copy assignment
        void operator=(const multiRegionCoupledDecomp&) = delete;


    //- Runtime type information
    TypeName("multiRegionCoupled");


    // Constructors

        //- Construct given decomposition dictionary and optional region name
        explicit multiRegionCoupledDecomp
        (
            const dictionary& decompDict,
            const word& regionName = ""
        );


    virtual ~multiRegionCoupledDecomp() = default;


    // Member Functions

        //- Serial only
        virtual bool parallelAware() const
        {
            return false;
        }

        //- Return the processor for every cell of the region mesh
        virtual labelList decompose
        (
            const polyMesh& mesh,
            const pointField& cc,
            const scalarField& cWeights = scalarField::null()
        ) const;

        //- Not implemented: the weld is built from the region meshes
        virtual labelList decompose
        (
            const CompactListList<label>& globalCellCells,
            const pointField& cc,
            const scalarField& cWeights = scalarField::null()
        ) const
        {
            NotImplemented;
            return labelList();
        }

        //- Not implemented: the weld is built from the region meshes
        virtual labelList decompose
        (
            const labelListList& globalCellCells,
            const pointField& cc,
            const scalarField& cWeights = scalarField::null()
        ) const
        {
            NotImplemented;
            return labelList();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
