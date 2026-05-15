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
    Foam::regionModels::areaSurfaceFilmModels::circular

Description
    Circular footprint shape.

    Builds a deterministic tangent-plane basis from the source-face normal,
    independent of any user-supplied direction, and sets equal major and
    minor footprint radii, i.e. an isotropic footprint.

Usage
    Minimal example:
    \verbatim
    physicalSpreading
    {
        // Optional entries
        footprintShape  circular;
    }
    \endverbatim

    where the entries mean:
    \table
      Property       | Description                     | Type | Reqd | Deflt
      footprintShape | Type name: circular             | word | yes  | -
    \endtable

SourceFiles
    circular.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_circular_H
#define regionFaModels_circular_H

#include "footprintShape.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                          Class circular Declaration
\*---------------------------------------------------------------------------*/

class circular
:
    public footprintShape
{
public:

    //- Runtime type information
    TypeName("circular");


    // Constructors

        //- Construct from dictionary
        explicit circular(const dictionary& dict);


    //- Destructor
    virtual ~circular() = default;


    // Member Functions

        //- Build an isotropic tangent-plane footprint for a source face
        virtual footprintGeometry footprintAt
        (
            const physicalSpreadingGeometry& geom,
            const label sourceGlobalFace,
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
