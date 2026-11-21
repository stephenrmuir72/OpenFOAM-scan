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
    Foam::regionModels::areaSurfaceFilmModels::elliptical

Description
    Elliptical footprint shape.

    Projects a user-supplied \c direction onto the source-face tangent
    plane to build the major axis, and scales the minor radius by
    \c aspectRatio to obtain the major radius, i.e. an anisotropic
    footprint aligned with \c direction.

Usage
    Minimal example:
    \verbatim
    spatialSpreading
    {
        // Optional entries
        footprintShape    elliptical;

        // Mandatory entries
        aspectRatio       <scalar>;    // Major/minor axis ratio, >= 1
        direction         <vector>;    // Major-axis direction
    }
    \endverbatim

    where the entries mean:
    \table
      Property    | Description                      | Type   | Reqd | Deflt
      footprintShape | Type name: elliptical         | word   | yes  | -
      aspectRatio | Ellipse major/minor ratio, >= 1  | scalar | yes  | -
      direction   | Ellipse major-axis direction     | vector | yes  | -
    \endtable

SourceFiles
    elliptical.cxx

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_elliptical_H
#define regionFaModels_elliptical_H

#include "footprintShape.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                         Class elliptical Declaration
\*---------------------------------------------------------------------------*/

class elliptical
:
    public footprintShape
{
    // Private Data

        //- Ellipse major-to-minor axis ratio
        const scalar aspectRatio_;

        //- Ellipse major-axis direction
        const vector direction_;


public:

    //- Runtime type information
    TypeName("elliptical");


    // Constructors

        //- Construct from dictionary
        explicit elliptical(const dictionary& dict);


    //- Destructor
    virtual ~elliptical() = default;


    // Member Functions

        //- Build an anisotropic tangent-plane footprint for a source face
        virtual footprintGeometry footprintAt
        (
            const spatialSpreadingData& data,
            const label sourceFacei,
            const scalar minorRadius
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
