/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2018 OpenCFD Ltd.
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

#include "multiBandAbsorption.H"
#include "addToRunTimeSelectionTable.H"
#include "faceZoneOverrides.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(multiBandAbsorption, 0);

        addToRunTimeSelectionTable
        (
            wallAbsorptionEmissionModel,
            multiBandAbsorption,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * //

namespace Foam
{
namespace
{

//- Per-band coefficients for each zone block that sets keyName, with the
//  patch-local index of the last zone block that sets it
void readZones
(
    const polyPatch& pp,
    const dictionary& zones,
    const Foam::radiation::faceZoneOverrides::selection& sel,
    const word& keyName,
    const label nBands,
    labelList& zoneIndex,
    List<scalarList>& zoneCoeffs
)
{
    if (!sel.setsAny(keyName))
    {
        return;
    }

    zoneIndex = sel.lastZone(keyName);

    const wordList& names = sel.names();

    zoneCoeffs.setSize(names.size());

    forAll(names, j)
    {
        if (sel.sets(j, keyName))
        {
            scalarList coeffs;
            zones.subDict(names[j]).readEntry(keyName, coeffs);

            if (coeffs.size() != nBands)
            {
                FatalErrorInFunction
                    << "faceZone '" << names[j] << "' sets '" << keyName
                    << "' with " << coeffs.size() << " bands, expected "
                    << nBands << " (the number of bands of the base list)"
                    << " on patch " << pp.name() << exit(FatalError);
            }

            zoneCoeffs[j] = coeffs;
        }
    }
}


//- Values of one band: the base value on every face, then the zone values
//  on the faces each zone selects
tmp<scalarField> bandValues
(
    const label size,
    const scalarList& base,
    const labelList& zoneIndex,
    const List<scalarList>& zoneCoeffs,
    const label bandI
)
{
    tmp<scalarField> tvalues = tmp<scalarField>::New(size, base[bandI]);

    if (!zoneIndex.empty())
    {
        scalarField& values = tvalues.ref();

        forAll(zoneIndex, i)
        {
            const label j = zoneIndex[i];

            if (j >= 0)
            {
                values[i] = zoneCoeffs[j][bandI];
            }
        }
    }

    return tvalues;
}

}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::multiBandAbsorption::multiBandAbsorption
(
    const dictionary& dict,
    const polyPatch& pp
)
:
    wallAbsorptionEmissionModel(dict, pp),
    coeffsDict_(dict),
    aCoeffs_(),
    eCoeffs_(),
    nBands_(0),
    nFaces_(pp.size()),
    aZone_(),
    aZoneCoeffs_(),
    eZone_(),
    eZoneCoeffs_()
{
    coeffsDict_.readEntry("absorptivity", aCoeffs_);
    coeffsDict_.readEntry("emissivity", eCoeffs_);

    if (eCoeffs_.size() != aCoeffs_.size())
    {
        FatalErrorInFunction
            << "Number of absorptivity bands (" << aCoeffs_.size()
            << ") must equal the number of emissivity bands ("
            << eCoeffs_.size() << ") on patch " << pp.name()
            << exit(FatalError);
    }

    nBands_ = aCoeffs_.size();

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

    faceZoneOverrides::checkKeys(zones, allowed, false);

    const faceZoneOverrides::selection sel(pp, zones, allowed);

    readZones(pp, zones, sel, "absorptivity", nBands_, aZone_, aZoneCoeffs_);
    readZones(pp, zones, sel, "emissivity", nBands_, eZone_, eZoneCoeffs_);

    sel.overlapWarning();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::multiBandAbsorption::~multiBandAbsorption()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::radiation::multiBandAbsorption::a
(
    const label bandI,
    const vectorField* incomingDirection,
    const scalarField* T
) const
{
    return bandValues(nFaces_, aCoeffs_, aZone_, aZoneCoeffs_, bandI);
}

Foam::scalar Foam::radiation::multiBandAbsorption::a
(
    const label faceI,
    const label bandI,
    const vector dir,
    const scalar T
) const
{
    const label j = aZone_.empty() ? -1 : aZone_[faceI];

    return j >= 0 ? aZoneCoeffs_[j][bandI] : aCoeffs_[bandI];
}


Foam::tmp<Foam::scalarField> Foam::radiation::multiBandAbsorption::e
(
    const label bandI,
    const vectorField* incomingDirection,
    const scalarField* T
) const
{
    return bandValues(nFaces_, eCoeffs_, eZone_, eZoneCoeffs_, bandI);
}


Foam::scalar Foam::radiation::multiBandAbsorption::e
(
    const label faceI,
    const label bandI,
    const vector dir,
    const scalar T
) const
{
    const label j = eZone_.empty() ? -1 : eZone_[faceI];

    return j >= 0 ? eZoneCoeffs_[j][bandI] : eCoeffs_[bandI];
}


// ************************************************************************* //
