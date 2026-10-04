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
    Foam::regionModels::areaSurfaceFilmModels::footprintShape

Description
    Abstract base class for the runtime-selectable footprint-shape
    submodel of the \c spatialSpreading model.

    A footprint-shape model builds the source-face tangent plane and the
    major/minor footprint radii for a scaled parcel radius, and provides the
    shared tangent-plane projection and ellipse-boundary-radius evaluation used
    by every option regardless of orientation.

Usage
    Minimal example:
    \verbatim
    spatialSpreading
    {
        // Optional entries
        footprintShape    <word>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property       | Description             | Type | Reqd  | Deflt
      footprintShape | Type name               | word | no    | circular
    \endtable

    Options for the \c footprintShape entry:
    \verbatim
      circular    | Use circular footprint shape (default)
      elliptical  | Use elliptical footprint shape
    \endverbatim

SourceFiles
    footprintShape.cxx
    footprintShapeNew.cxx

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_footprintShape_H
#define regionFaModels_footprintShape_H

#include "autoPtr.H"
#include "dictionary.H"
#include "point.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// Forward Declarations
class spatialSpreadingData;

/*---------------------------------------------------------------------------*\
                      Class footprintShape Declaration
\*---------------------------------------------------------------------------*/

class footprintShape
{
protected:

    // Protected Data

        //- Const reference to model dictionary
        const dictionary& dict_;


public:

    // Public Data Types

        //- Compact tangent-plane description of a source footprint
        struct footprintGeometry
        {
            //- Source face centre
            point centre = Foam::zero{};

            //- Source face unit normal
            vector normal = Foam::zero{};

            //- Unit major axis within the source tangent plane
            vector majorAxis = Foam::zero{};

            //- Unit minor axis within the source tangent plane
            vector minorAxis = Foam::zero{};

            //- Footprint major-axis radius
            scalar majorRadius = 0;

            //- Footprint minor-axis radius
            scalar minorRadius = 0;

            //- Conservative bounding radius, max(majorRadius, minorRadius)
            scalar enclosingRadius() const
            {
                return Foam::max(majorRadius, minorRadius);
            }
        };


    //- Runtime type information
    TypeName("footprintShape");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            footprintShape,
            dictionary,
            (
                const dictionary& dict
            ),
            (dict)
        );


    // Constructors

        //- Construct from model type and dictionary
        footprintShape(const word& modelType, const dictionary& dict);


    // Generated Methods

        //- No copy construct
        footprintShape(const footprintShape&) = delete;

        //- No copy assignment
        void operator=(const footprintShape&) = delete;


    // Selectors

        //- Select a footprint-shape model
        static autoPtr<footprintShape> New(const dictionary& dict);


    //- Destructor
    virtual ~footprintShape() = default;


    // Member Functions

        //- Build the tangent-plane footprint geometry for a source face
        //- given the scaled (beta*d) minor radius
        virtual footprintGeometry footprintAt
        (
            const spatialSpreadingData& data,
            const label sourceFacei,
            const scalar minorRadius
        ) const = 0;

        //- Tangent-plane (x,y) offset of a point relative to a footprint
        static void planarOffset
        (
            const footprintGeometry& fp,
            const point& target,
            scalar& x,
            scalar& y
        );

        //- Direction-dependent ellipse/circle boundary radius at tangent-plane
        scalar footprintRadius
        (
            const scalar x,
            const scalar y,
            const scalar majorRadius,
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
