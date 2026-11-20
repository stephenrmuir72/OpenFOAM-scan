/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2019 OpenCFD Ltd.
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

#include "lookup.H"
#include "addToRunTimeSelectionTable.H"
#include "faceZoneOverrides.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(lookup, 0);
        addToRunTimeSelectionTable
        (
            boundaryRadiationPropertiesPatch,
            lookup,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * //

namespace Foam
{
namespace
{

//- Scalar value of keyName. Dictionary values are rejected
scalar readScalar
(
    const dictionary& dict,
    const word& keyName,
    const polyPatch& pp
)
{
    if (dict.isDict(keyName))
    {
        FatalErrorInFunction
            << "'" << keyName << "' on patch " << pp.name()
            << " must be a scalar. lookup does not accept a dictionary "
            << "(e.g. mappedFile) value" << exit(FatalError);
    }

    return dict.get<scalar>(keyName);
}


//- Per-face values of keyName: the base value, overwritten by the zone value
//  of the last zone that selects each face. Left empty if no zone sets
//  keyName, so the base value is used everywhere. faceSet is true for each
//  face that has a value
void readFaceValues
(
    const polyPatch& pp,
    const dictionary& zones,
    const Foam::radiation::faceZoneOverrides::selection& sel,
    const word& keyName,
    const scalar base,
    const bool baseSet,
    scalarField& faceValues,
    boolList& faceSet
)
{
    if (!sel.setsAny(keyName))
    {
        return;
    }

    const labelList last = sel.lastZone(keyName);
    const wordList& names = sel.names();

    scalarList zoneValues(names.size(), 0);

    forAll(names, j)
    {
        if (sel.sets(j, keyName))
        {
            zoneValues[j] =
                readScalar(zones.subDict(names[j]), keyName, pp);
        }
    }

    faceValues = scalarField(pp.size(), base);
    faceSet = boolList(pp.size(), baseSet);

    forAll(last, i)
    {
        if (last[i] >= 0)
        {
            faceValues[i] = zoneValues[last[i]];
            faceSet[i] = true;
        }
    }
}

}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::lookup::lookup
(
    const dictionary& dict,
    const polyPatch& pp
)
:
    boundaryRadiationPropertiesPatch(dict, pp),
    pp_(pp),
    e0_(readScalar(dict, "emissivity", pp)),
    a0_(dict.found("absorptivity") ? readScalar(dict, "absorptivity", pp) : 0),
    t0_
    (
        dict.found("transmissivity")
      ? readScalar(dict, "transmissivity", pp)
      : 0
    ),
    aSet_(dict.found("absorptivity")),
    nAMissing_(0),
    eFace_(),
    aFace_(),
    tFace_(),
    aFaceSet_()
{
    const dictionary noZones;

    const dictionary& zones =
        dict.found("faceZoneOverrides")
      ? dict.subDict("faceZoneOverrides")
      : noZones;

    const wordList allowed
    {
        "emissivity",
        "absorptivity",
        "transmissivity"
    };

    faceZoneOverrides::checkKeys(zones, allowed, false);

    const faceZoneOverrides::selection sel(pp, zones, allowed);

    boolList eSet;
    boolList tSet;

    readFaceValues(pp, zones, sel, "emissivity", e0_, true, eFace_, eSet);
    readFaceValues
    (
        pp, zones, sel, "absorptivity", a0_, aSet_, aFace_, aFaceSet_
    );
    readFaceValues(pp, zones, sel, "transmissivity", t0_, true, tFace_, tSet);

    if (aFace_.empty())
    {
        nAMissing_ = aSet_ ? 0 : pp.size();
    }
    else
    {
        label nASet = 0;

        forAll(aFaceSet_, faceI)
        {
            nASet += aFaceSet_[faceI];
        }

        nAMissing_ = pp.size() - nASet;
    }

    sel.overlapWarning();
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::radiation::lookup::e
(
    const label bandI,
    const vectorField* dir,
    const scalarField* T
) const
{
    if (eFace_.empty())
    {
        return tmp<scalarField>::New(pp_.size(), e0_);
    }

    return tmp<scalarField>::New(eFace_);
}


Foam::scalar Foam::radiation::lookup::e
(
    const label faceI,
    const label bandI,
    const vector& dir,
    const scalar T
) const
{
    return eFace_.empty() ? e0_ : eFace_[faceI];
}


Foam::tmp<Foam::scalarField> Foam::radiation::lookup::a
(
    const label bandI,
    const vectorField* dir,
    const scalarField* T
) const
{
    if (nAMissing_)
    {
        FatalErrorInFunction
            << "absorptivity is not set for " << nAMissing_ << " face(s) on "
            << "patch " << pp_.name() << ". Set absorptivity in the patch "
            << "entry or in a faceZoneOverrides block that selects those "
            << "faces" << exit(FatalError);
    }

    if (aFace_.empty())
    {
        return tmp<scalarField>::New(pp_.size(), a0_);
    }

    return tmp<scalarField>::New(aFace_);
}


Foam::scalar Foam::radiation::lookup::a
(
    const label faceI,
    const label bandI,
    const vector& dir,
    const scalar T
) const
{
    if (aFace_.empty())
    {
        if (!aSet_)
        {
            FatalErrorInFunction
                << "absorptivity is not set on patch " << pp_.name()
                << exit(FatalError);
        }

        return a0_;
    }

    if (!aFaceSet_[faceI])
    {
        FatalErrorInFunction
            << "absorptivity is not set for face " << faceI << " on patch "
            << pp_.name() << exit(FatalError);
    }

    return aFace_[faceI];
}


Foam::tmp<Foam::scalarField> Foam::radiation::lookup::t
(
    const label bandI,
    const vectorField* dir,
    const scalarField* T
) const
{
    if (tFace_.empty())
    {
        return tmp<scalarField>::New(pp_.size(), t0_);
    }

    return tmp<scalarField>::New(tFace_);
}


Foam::scalar Foam::radiation::lookup::t
(
    const label faceI,
    const label bandI,
    const vector& dir,
    const scalar T
) const
{
    return tFace_.empty() ? t0_ : tFace_[faceI];
}


Foam::tmp<Foam::scalarField> Foam::radiation::lookup::rSpec
(
    const label bandI,
    const vectorField* dir,
    const scalarField* T
) const
{
    return tmp<scalarField>::New(pp_.size(), Zero);
}


Foam::scalar Foam::radiation::lookup::rSpec
(
    const label faceI,
    const label bandI,
    const vector& dir,
    const scalar T
) const
{
    return Zero;
}


Foam::tmp<Foam::scalarField> Foam::radiation::lookup::rDiff
(
    const label bandI,
    const vectorField* dir,
    const scalarField* T
) const
{
    return tmp<scalarField>::New(pp_.size(), Zero);
}


Foam::scalar Foam::radiation::lookup::rDiff
(
    const label faceI,
    const label bandI,
    const vector& dir,
    const scalar T
) const
{
    return Zero;
}


bool Foam::radiation::lookup::isGrey() const
{
    return true;
}


Foam::label Foam::radiation::lookup::nBands() const
{
    return 1;
}


// ************************************************************************* //
