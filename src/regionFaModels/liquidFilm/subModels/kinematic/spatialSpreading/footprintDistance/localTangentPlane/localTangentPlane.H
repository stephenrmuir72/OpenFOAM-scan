/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Keysight Technologies
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
    Foam::regionModels::areaSurfaceFilmModels::localTangentPlane

Description
    Local tangent-plane footprint distance.

    Finds candidates using the local face-centre search tree and
    outputs the straight-line tangent-plane distance
    \c sqrt(x^2 + y^2) as the selected-distance metric for every candidate.

Usage
    Minimal example:
    \verbatim
    spatialSpreading
    {
        // Optional entries
        footprintDistance  localTangentPlane;
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Reqd | Deflt
      footprintDistance | Type name: localTangentPlane  | word | yes  | -
    \endtable

SourceFiles
    localTangentPlane.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_localTangentPlane_H
#define regionFaModels_localTangentPlane_H

#include "footprintDistance.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                     Class localTangentPlane Declaration
\*---------------------------------------------------------------------------*/

class localTangentPlane
:
    public footprintDistance
{
public:

    //- Runtime type information
    TypeName("localTangentPlane");


    // Constructors

        //- Construct from dictionary
        explicit localTangentPlane(const dictionary& dict);


    //- Destructor
    virtual ~localTangentPlane() = default;


    // Member Functions

        //- Find candidates using the local search tree
        virtual void findCandidates
        (
            const spatialSpreadingData& data,
            const footprintShape& shape,
            const footprintShape::footprintGeometry& fp,
            const label sourceProc,
            const label sourceFace,
            const scalar searchRadius,
            DynamicList<label>& candidateFaces,
            DynamicList<scalar>& selectedDistances
        ) const;

        //- Return the projected-centroid tangent-plane distance unchanged
        virtual scalar kernelDistance
        (
            const scalar selectedDistance,
            const scalar planarDistance,
            const scalar centroidPlanarDistance
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
