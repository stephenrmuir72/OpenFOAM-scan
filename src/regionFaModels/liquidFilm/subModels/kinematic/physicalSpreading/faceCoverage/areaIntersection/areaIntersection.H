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
    Foam::regionModels::areaSurfaceFilmModels::areaIntersection

Description
    Polygon-intersection coverage measure.

    Projects the target face polygon onto the footprint tangent plane and
    computes the exact overlap area and centroid with the footprint polygon
    using \c faceAreaIntersect, giving partial-face weights. The returned
    coverage result contains the overlap area and its centroid.

    Requires the shared global face polygons to be assembled.

Usage
    Minimal example:
    \verbatim
    physicalSpreading
    {
        // Optional entries
        faceCoverage  areaIntersection;
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Reqd | Deflt
      faceCoverage | Type name: areaIntersection        | word | yes  | -
    \endtable

SourceFiles
    areaIntersection.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_areaIntersection_H
#define regionFaModels_areaIntersection_H

#include "faceCoverage.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                       Class areaIntersection Declaration
\*---------------------------------------------------------------------------*/

class areaIntersection
:
    public faceCoverage
{
    // Private Member Functions

        //- Polygonal footprint/face overlap area and centroid in the
        //- source tangent plane
        scalar overlapArea
        (
            const physicalSpreadingGeometry& geom,
            const footprintShape::footprintGeometry& fp,
            const label targetGlobalFace,
            point& overlapCentre
        ) const;


public:

    //- Runtime type information
    TypeName("areaIntersection");


    // Constructors

        //- Construct from dictionary
        explicit areaIntersection(const dictionary& dict);


    //- Destructor
    virtual ~areaIntersection() = default;


    // Member Functions

        //- Requires the shared global face polygons to be assembled
        virtual bool needsFacePoints() const noexcept { return true; }

        //- The overlap centroid may differ from the target face centre
        virtual bool usesProjectedCentroid() const noexcept { return true; }

        //- Conservative broad-phase reject, then exact polygon overlap
        virtual coverageResult evaluate
        (
            const physicalSpreadingGeometry& geom,
            const footprintShape::footprintGeometry& fp,
            const label targetGlobalFace,
            const scalar selectedDistance,
            const scalar boundaryRadius
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
