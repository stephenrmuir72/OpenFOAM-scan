/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
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
    Foam::cellVolumeWeightMethod

Description
    Cell-volume-weighted mesh-to-mesh interpolation class

    Piecewise-constant volume transfer with complete box-candidate searches.
    Small overlaps and pair bounds use the smaller participating geometric
    volume: the same oriented first-vertex fan as the intersection engine.
    This can differ from FV volume on warped faces. Volume-carrying fan tets
    must be positive; exactly collinear face triangles are ignored. Numerically
    unresolved tetrahedron areas are errors.

    Conservative on complete, non-overlapping partitions of the represented
    geometry, up to discarded overlaps and geometric roundoff. This is not an
    FV-volume conservation guarantee when the representations differ.
    Normalization of partially covered cells changes that guarantee.

SourceFiles
    cellVolumeWeightMethod.C

\*---------------------------------------------------------------------------*/

#ifndef cellVolumeWeightMethod_H
#define cellVolumeWeightMethod_H

#include "meshToMeshMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class cellVolumeWeightMethod Declaration
\*---------------------------------------------------------------------------*/

class cellVolumeWeightMethod
:
    public meshToMeshMethod
{
protected:

    // Protected Data
        label candidateTests_;
        scalar discardedVolume_;

    // Protected Member Functions
        //- Calculate the mesh-to-mesh addressing and weights
        void calculateAddressing
        (
            labelListList& srcToTgtCellAddr,
            scalarListList& srcToTgtCellWght,
            pointListList& srcToTgtCellVec,
            labelListList& tgtToSrcCellAddr,
            scalarListList& tgtToSrcCellWght,
            pointListList& tgtToSrcCellVec,
            const bool centroids
        );


private:

    // Private member functions

        //- No copy construct
        cellVolumeWeightMethod(const cellVolumeWeightMethod&) = delete;

        //- No copy assignment
        void operator=(const cellVolumeWeightMethod&) = delete;


public:

    //- Run-time type information
    TypeName("cellVolumeWeight");

    //- Construct from source and target meshes
    cellVolumeWeightMethod(const polyMesh& src, const polyMesh& tgt);

    //- Destructor
    virtual ~cellVolumeWeightMethod();


    // Member Functions

        label candidateTests() const noexcept { return candidateTests_; }
        scalar discardedVolume() const noexcept { return discardedVolume_; }

        // Evaluate

            //- Calculate addressing and weights and optionally offset vectors
            virtual void calculate
            (
                labelListList& srcToTgtAddr,
                scalarListList& srcToTgtWght,
                pointListList& srcToTgtVec,
                labelListList& tgtToSrcAddr,
                scalarListList& tgtToSrcWght,
                pointListList& tgtToSrcVec
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
