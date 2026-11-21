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

\*---------------------------------------------------------------------------*/

#include "elliptical.H"
#include "addToRunTimeSelectionTable.H"
#include "spatialSpreadingData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(elliptical, 0);
addToRunTimeSelectionTable(footprintShape, elliptical, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

elliptical::elliptical(const dictionary& dict)
:
    footprintShape(typeName, dict),
    aspectRatio_(dict.get<scalar>("aspectRatio")),
    direction_(dict.get<vector>("direction"))
{
    if (!std::isfinite(aspectRatio_) || !(aspectRatio_ >= 1))
    {
        FatalIOErrorInFunction(dict)
            << "aspectRatio must be finite and >= 1 for an elliptical "
            << "footprint, received " << aspectRatio_
            << exit(FatalIOError);
    }

    const scalar directionMagnitude = Foam::mag(direction_);
    if
    (
        !std::isfinite(directionMagnitude)
     || directionMagnitude <= ROOTVSMALL
    )
    {
        FatalIOErrorInFunction(dict)
            << "direction must be finite and non-zero for an elliptical "
            << "footprint"
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

footprintShape::footprintGeometry elliptical::footprintAt
(
    const spatialSpreadingData& data,
    const label sourceFacei,
    const scalar minorRadius
) const
{
    footprintGeometry fp;
    fp.centre = data.centres()[sourceFacei];
    fp.normal = data.normals()[sourceFacei];

    if (aspectRatio_ == 1)
    {
        vector reference(1, 0, 0);
        if (Foam::mag(fp.normal.y()) < Foam::mag(fp.normal.x()))
        {
            reference = vector(0, 1, 0);
        }
        if
        (
            Foam::mag(fp.normal.z()) < Foam::mag(fp.normal.x())
         && Foam::mag(fp.normal.z()) < Foam::mag(fp.normal.y())
        )
        {
            reference = vector(0, 0, 1);
        }

        fp.majorAxis = reference;
        fp.majorAxis.removeCollinear(fp.normal);
    }
    else
    {
        fp.majorAxis = direction_;
        fp.majorAxis.removeCollinear(fp.normal);
    }

    const scalar majorAxisMagnitude = Foam::mag(fp.majorAxis);
    if
    (
        !std::isfinite(majorAxisMagnitude)
     || majorAxisMagnitude <= ROOTVSMALL
    )
    {
        if (aspectRatio_ == 1)
        {
            FatalIOErrorInFunction(dict_)
                << "Cannot construct a finite tangent-plane basis for "
                << "source face " << sourceFacei
                << exit(FatalIOError);
        }
        else
        {
            FatalIOErrorInFunction(dict_)
                << "direction " << direction_
                << " is parallel to the normal of active source face "
                << sourceFacei
                << "; its tangent-plane projection is degenerate"
                << exit(FatalIOError);
        }
    }

    fp.majorAxis.normalise();
    fp.minorAxis = fp.normal.cross(fp.majorAxis);
    const scalar minorAxisMagnitude = Foam::mag(fp.minorAxis);

    if (!std::isfinite(minorAxisMagnitude) || minorAxisMagnitude <= ROOTVSMALL)
    {
        FatalIOErrorInFunction(dict_)
            << "Cannot construct a finite tangent-plane basis for source "
            << "face " << sourceFacei
            << exit(FatalIOError);
    }
    fp.minorAxis.normalise();

    fp.majorRadius = aspectRatio_*minorRadius;
    if (!std::isfinite(fp.majorRadius))
    {
        FatalErrorInFunction
            << "Non-finite major radius for aspectRatio " << aspectRatio_
            << " and minor radius " << minorRadius
            << " on source face " << sourceFacei
            << abort(FatalError);
    }

    fp.minorRadius = minorRadius;

    return fp;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
