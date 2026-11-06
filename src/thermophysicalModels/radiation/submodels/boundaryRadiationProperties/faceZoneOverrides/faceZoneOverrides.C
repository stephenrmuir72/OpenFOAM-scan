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

#include "faceZoneOverrides.H"
#include "polyMesh.H"
#include "faceZone.H"
#include "processorPolyPatch.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * Functions * * * * * * * * * * * * * * * //

bool Foam::radiation::faceZoneOverrides::present(const dictionary& dict)
{
    if (dict.found("faceZoneOverrides"))
    {
        return true;
    }

    return
        dict.isDict("wallAbsorptionEmissionModel")
     && dict.subDict("wallAbsorptionEmissionModel").found("faceZoneOverrides");
}


void Foam::radiation::faceZoneOverrides::checkKeys
(
    const dictionary& overrides,
    const wordList& allowed,
    const bool allowCoeffs
)
{
    for (const word& zoneName : overrides.toc())
    {
        if (!overrides.isDict(zoneName))
        {
            FatalErrorInFunction
                << "faceZoneOverrides entry '" << zoneName
                << "' must be a dictionary" << exit(FatalError);
        }

        const dictionary& zone = overrides.subDict(zoneName);

        for (const word& key : zone.toc())
        {
            if (allowed.found(key))
            {
                continue;
            }

            const word coeffsSuffix("Coeffs");

            const word base =
                key.size() > coeffsSuffix.size()
             && key.substr(key.size() - coeffsSuffix.size()) == coeffsSuffix
              ? word(key.substr(0, key.size() - coeffsSuffix.size()))
              : word::null;

            if (allowCoeffs && allowed.found(base))
            {
                if (!zone.found(base))
                {
                    FatalErrorInFunction
                        << "Key '" << key << "' in faceZoneOverrides block '"
                        << zoneName << "' needs '" << base << "' in the same "
                        << "block" << exit(FatalError);
                }

                continue;
            }

            FatalErrorInFunction
                << "Unknown key '" << key << "' in faceZoneOverrides "
                << "block '" << zoneName << "'" << nl
                << "Allowed keys: " << allowed << exit(FatalError);
        }
    }
}


Foam::radiation::faceZoneOverrides::selection::selection
(
    const polyPatch& pp,
    const dictionary& overrides,
    const wordList& keys
)
:
    patchName_(pp.name()),
    nFaces_(pp.size()),
    isProcessor_(isA<processorPolyPatch>(pp)),
    keys_(keys)
{
    const faceZoneMesh& fzs = pp.boundaryMesh().mesh().faceZones();

    const label patchStart = pp.start();


    DynamicList<word> names;
    DynamicList<wordList> zoneKeys;
    DynamicList<labelList> faces;

    for (const word& zoneName : overrides.toc())
    {
        const dictionary& zone = overrides.subDict(zoneName);

        DynamicList<word> setKeys;

        for (const word& key : keys)
        {
            if (zone.found(key))
            {
                setKeys.append(key);
            }
        }

        if (setKeys.empty())
        {
            continue;
        }

        const label zonei = fzs.findZoneID(zoneName);

        if (zonei < 0)
        {
            FatalErrorInFunction
                << "Unknown faceZone '" << zoneName
                << "' in faceZoneOverrides on patch " << pp.name()
                << exit(FatalError);
        }

        const faceZone& fz = fzs[zonei];

        DynamicList<label> localFaces;

        forAll(fz, i)
        {
            const label facei = fz[i];

            if (facei >= patchStart && facei < patchStart + nFaces_)
            {
                localFaces.append(facei - patchStart);
            }
        }

        labelList selected(localFaces);
        inplaceUniqueSort(selected);

        if (!isProcessor_)
        {
            const label nGlobal = returnReduce(selected.size(), sumOp<label>());

            if (!nGlobal)
            {
                WarningInFunction
                    << "faceZone '" << zoneName << "' selects no faces on "
                    << "patch " << pp.name() << " for " << wordList(setKeys)
                    << endl;
            }
        }

        names.append(zoneName);
        zoneKeys.append(wordList(setKeys));
        faces.append(std::move(selected));
    }

    names_ = wordList(names);
    zoneKeys_ = List<wordList>(zoneKeys);
    faces_ = List<labelList>(faces);
}


bool Foam::radiation::faceZoneOverrides::selection::sets
(
    const label j,
    const word& keyName
) const
{
    return zoneKeys_[j].found(keyName);
}


bool Foam::radiation::faceZoneOverrides::selection::setsAny
(
    const word& keyName
) const
{
    forAll(names_, j)
    {
        if (sets(j, keyName))
        {
            return true;
        }
    }

    return false;
}


Foam::labelList Foam::radiation::faceZoneOverrides::selection::lastZone
(
    const word& keyName
) const
{
    labelList last(nFaces_, -1);

    forAll(names_, j)
    {
        if (sets(j, keyName))
        {
            for (const label i : faces_[j])
            {
                last[i] = j;
            }
        }
    }

    return last;
}


void Foam::radiation::faceZoneOverrides::selection::overlapWarning() const
{
    if (isProcessor_ || names_.empty())
    {
        return;
    }

    labelList shared(nFaces_, 0);

    forAll(keys_, k)
    {
        labelList nZones(nFaces_, 0);

        forAll(names_, j)
        {
            if (sets(j, keys_[k]))
            {
                for (const label i : faces_[j])
                {
                    ++nZones[i];
                }
            }
        }

        forAll(nZones, i)
        {
            if (nZones[i] > 1)
            {
                shared[i] = 1;
            }
        }
    }

    label nShared = 0;

    forAll(shared, i)
    {
        nShared += shared[i];
    }

    nShared = returnReduce(nShared, sumOp<label>());

    if (nShared)
    {
        WarningInFunction
            << nShared << " faces on patch " << patchName_ << " are selected "
            << "by more than one faceZone block for the same property." << nl
            << "The last faceZone block in dictionary order is used." << endl;
    }
}

// ************************************************************************* //
