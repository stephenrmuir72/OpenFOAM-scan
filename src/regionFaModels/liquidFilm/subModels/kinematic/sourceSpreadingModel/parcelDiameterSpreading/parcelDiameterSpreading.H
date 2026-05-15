/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 OpenCFD Ltd.
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
    Foam::regionModels::areaSurfaceFilmModels::parcelDiameterSpreading

Description
    Spreads Lagrangian-to-film source terms over a circular footprint
    whose diameter equals the impacting Lagrangian parcel diameter.

    For each finite-area face carrying a non-zero source, the spreading
    radius is \c diameter/2.  Neighbouring faces whose centroids lie within
    this radius are found by breadth-first search (BFS) through
    \c PrimitivePatch::faceFaces() adjacency.  Source is distributed
    proportionally to face area (conservative).

Usage
    \verbatim
    sourceSpreading
    {
        type    parcelDiameter;
    }
    \endverbatim

Note
    - The parcel diameter must be provided via
      \c liquidFilmBase::addDiameterSource() from the Lagrangian side.
    - On processor boundaries, only processor-local faces are included
      in the stencil.

SourceFiles
    parcelDiameterSpreading.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_parcelDiameterSpreading_H
#define regionFaModels_parcelDiameterSpreading_H

#include "sourceSpreadingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                  Class parcelDiameterSpreading Declaration
\*---------------------------------------------------------------------------*/

class parcelDiameterSpreading
:
    public sourceSpreadingModel
{
    // Private Member Functions

        //- Templated implementation of spread
        template<class Type>
        void spreadImpl
        (
            Field<Type>& field,
            const scalarField& diameter
        ) const;

        //- No copy construct
        parcelDiameterSpreading(const parcelDiameterSpreading&) = delete;

        //- No copy assignment
        void operator=(const parcelDiameterSpreading&) = delete;


public:

    //- Runtime type information
    TypeName("parcelDiameter");


    // Constructors

        //- Construct from mesh and dictionary
        parcelDiameterSpreading(const faMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~parcelDiameterSpreading() = default;


    // Member Functions

        //- Spread scalar source field using parcel diameter footprint
        virtual void spread
        (
            scalarField& field,
            const scalarField& diameter
        ) const;

        //- Spread vector source field using parcel diameter footprint
        virtual void spread
        (
            vectorField& field,
            const scalarField& diameter
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
