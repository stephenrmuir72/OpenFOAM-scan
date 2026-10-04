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
    Foam::regionModels::areaSurfaceFilmModels::faceCentre

Description
    Face-centre coverage measure.

    Includes a candidate face when its selected distance lies within the
    footprint boundary radius, and uses the whole target face area as the
    coverage measure.

Usage
    Minimal example:
    \verbatim
    spatialSpreading
    {
        // Optional entries
        faceCoverage    faceCentre;
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Reqd | Deflt
      faceCoverage | Type name: faceCentre              | word | yes  | -
    \endtable

SourceFiles
    faceCentre.cxx

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_faceCentre_H
#define regionFaModels_faceCentre_H

#include "faceCoverage.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                         Class faceCentre Declaration
\*---------------------------------------------------------------------------*/

class faceCentre
:
    public faceCoverage
{
public:

    //- Runtime type information
    TypeName("faceCentre");


    // Constructors

        //- Construct from dictionary
        explicit faceCentre(const dictionary& dict);


    //- Destructor
    virtual ~faceCentre() = default;


    // Member Functions

        //- Calculate face-coverage data for a candidate face
        virtual coverageResult evaluate
        (
            const spatialSpreadingData& data,
            const footprintShape::footprintGeometry& fp,
            const label targetFacei,
            const scalar candidateDistance,
            const scalar footprintRadius
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
