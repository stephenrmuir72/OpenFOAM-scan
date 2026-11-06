/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2018, 2024 OpenCFD Ltd.
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

#include "constantAbsorption.H"
#include "addToRunTimeSelectionTable.H"
#include "faceZoneOverrides.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(constantAbsorption, 0);

        addToRunTimeSelectionTable
        (
            wallAbsorptionEmissionModel,
            constantAbsorption,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * //

namespace Foam
{
namespace
{

//- Zone functions for keyName, with the patch-face index of the last zone.
//  Left empty if no zone sets keyName
void readZones
(
    const polyPatch& pp,
    const dictionary& zones,
    const Foam::radiation::faceZoneOverrides::selection& sel,
    const word& keyName,
    labelList& zoneIndex,
    PtrList<Function1<scalar>>& zoneFuncs
)
{
    if (!sel.setsAny(keyName))
    {
        return;
    }

    zoneIndex = sel.lastZone(keyName);

    const wordList& names = sel.names();

    zoneFuncs.setSize(names.size());

    forAll(names, j)
    {
        if (sel.sets(j, keyName))
        {
            zoneFuncs.set
            (
                j,
                Function1<scalar>::New(keyName, zones.subDict(names[j]))
            );
        }
    }
}


//- Patch values from the base function on the faces no zone selects, then
//  the faceZone overrides. The base is not evaluated on overridden faces
tmp<scalarField> patchValues
(
    const Function1<scalar>& base,
    const label size,
    const scalarField* T,
    const labelList& zoneIndex,
    const PtrList<Function1<scalar>>& zoneFuncs,
    const word& keyName
)
{
    if (!size)
    {
        return tmp<scalarField>::New();
    }

    label nOverridden = 0;

    forAll(zoneIndex, i)
    {
        if (zoneIndex[i] >= 0)
        {
            ++nOverridden;
        }
    }

    // Without overrides, evaluate the base over the whole patch
    if (!nOverridden)
    {
        if (base.constant())
        {
            return tmp<scalarField>::New(size, base.value(0));
        }

        if (T)
        {
            return base.value(*T);
        }

        FatalErrorInFunction
            << "Attempted to set '" << keyName << "' using a "
            << "non-uniform function of Temperature, but temperature "
            << "field is unavailable" << abort(FatalError);

        return nullptr;
    }

    auto tvalues = tmp<scalarField>::New(size);
    scalarField& values = tvalues.ref();

    // The base is not needed when every face is overridden
    if (nOverridden < size)
    {
        if (base.constant())
        {
            const scalar baseValue = base.value(0);

            forAll(values, i)
            {
                if (zoneIndex[i] < 0)
                {
                    values[i] = baseValue;
                }
            }
        }
        else
        {
            if (!T)
            {
                FatalErrorInFunction
                    << "Attempted to set '" << keyName << "' using a "
                    << "non-uniform function of Temperature, but temperature "
                    << "field is unavailable" << abort(FatalError);
            }

            forAll(values, i)
            {
                if (zoneIndex[i] < 0)
                {
                    values[i] = base.value((*T)[i]);
                }
            }
        }
    }

    forAll(zoneIndex, i)
    {
        const label j = zoneIndex[i];

        if (j >= 0)
        {
            const Function1<scalar>& zoneFunc = zoneFuncs[j];

            if (zoneFunc.constant())
            {
                values[i] = zoneFunc.value(0);
            }
            else if (T)
            {
                values[i] = zoneFunc.value((*T)[i]);
            }
            else
            {
                FatalErrorInFunction
                    << "Attempted to set '" << keyName << "' in faceZone "
                    << "override using a non-uniform function of Temperature, "
                    << "but temperature field is unavailable"
                    << abort(FatalError);
            }
        }
    }

    return tvalues;
}


//- Single face value, from the faceZone override if it selects the face
scalar faceValue
(
    const Function1<scalar>& base,
    const labelList& zoneIndex,
    const PtrList<Function1<scalar>>& zoneFuncs,
    const label faceI,
    const scalar T
)
{
    if (zoneFuncs.size() && zoneIndex[faceI] >= 0)
    {
        return zoneFuncs[zoneIndex[faceI]].value(T);
    }

    return base.value(T);
}

}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::constantAbsorption::constantAbsorption
(
    const dictionary& dict,
     const polyPatch& pp
)
:
    wallAbsorptionEmissionModel(dict, pp),
    coeffsDict_(dict),
    a_(Function1<scalar>::New("absorptivity", coeffsDict_)),
    e_(Function1<scalar>::New("emissivity", coeffsDict_)),
    aZone_(),
    aZones_(),
    eZone_(),
    eZones_()
{
    const dictionary noZones;

    const dictionary& zones =
        dict.found("faceZoneOverrides")
      ? dict.subDict("faceZoneOverrides")
      : noZones;

    const wordList allowed
    {
        "absorptivity",
        "emissivity"
    };

    faceZoneOverrides::checkKeys(zones, allowed, true);

    const faceZoneOverrides::selection sel(pp, zones, allowed);

    readZones(pp, zones, sel, "absorptivity", aZone_, aZones_);
    readZones(pp, zones, sel, "emissivity", eZone_, eZones_);

    sel.overlapWarning();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::radiation::constantAbsorption::a
(
    const label bandI,
    const vectorField* incomingDirection,
    const scalarField* T
) const
{
    return patchValues
    (
        *a_,
        pp_.size(),
        T,
        aZone_,
        aZones_,
        "absorptivity"
    );
}


Foam::scalar Foam::radiation::constantAbsorption::a
(
    const label faceI,
    const label bandI,
    const vector dir,
    const scalar T
) const
{
    return faceValue(*a_, aZone_, aZones_, faceI, T);
}


Foam::tmp<Foam::scalarField> Foam::radiation::constantAbsorption::e
(
    const label bandI,
    const vectorField* incomingDirection,
    const scalarField* T
) const
{
    return patchValues
    (
        *e_,
        pp_.size(),
        T,
        eZone_,
        eZones_,
        "emissivity"
    );
}


Foam::scalar Foam::radiation::constantAbsorption::e
(
    const label faceI,
    const label bandI,
    const vector dir,
    const scalar T
) const
{
    return faceValue(*e_, eZone_, eZones_, faceI, T);
}


// ************************************************************************* //
