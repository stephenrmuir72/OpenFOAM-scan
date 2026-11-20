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

Namespace
    Foam::radiation::faceZoneOverrides

Description
    Helper functions for the optional faceZoneOverrides entry of patch
    radiation properties.

    faceZoneOverrides is a dictionary of faceZone blocks. Each block may set
    a subset of the property keys of the model. For a patch face, a property
    value is taken from the last zone block (in dictionary order) that
    selects the face and sets that property. Faces not selected by such a
    zone keep the patch value.

SourceFiles
    faceZoneOverrides.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_radiation_faceZoneOverrides_H
#define Foam_radiation_faceZoneOverrides_H

#include "dictionary.H"
#include "wordList.H"
#include "labelList.H"
#include "polyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{
namespace faceZoneOverrides
{

//- True if the entry, or its wallAbsorptionEmissionModel sub-dictionary,
//  contains a faceZoneOverrides entry
bool present(const dictionary& dict);

//- Fatal if any zone block of the overrides contains a key not in allowed.
//  If allowCoeffs, a '<key>Coeffs' entry is also accepted, provided <key> is
//  in allowed and is set in the same zone block
void checkKeys
(
    const dictionary& overrides,
    const wordList& allowed,
    const bool allowCoeffs
);


//- The zone blocks of a patch's faceZoneOverrides that set one of a set of
//  keys, with the patch-local faces of each zone block. Built once per patch
//  and shared by all keys. Warnings use counts summed over all ranks, so the
//  constructor and overlapWarning are collective on physical patches only.
//  Processor patches exist on some ranks only and skip the collectives
class selection
{
    // Private data

        //- Patch name, for messages
        word patchName_;

        //- Number of patch faces
        label nFaces_;

        //- True for a processor patch, which is not a collective participant
        bool isProcessor_;

        //- Keys that the selection covers
        wordList keys_;

        //- Zone block names that set a key, in dictionary order
        wordList names_;

        //- Keys of the covered set that each zone block sets
        List<wordList> zoneKeys_;

        //- Patch-local faces of each zone block: sorted, without duplicates
        List<labelList> faces_;


public:

    //- Construct from the zone blocks of overrides that set a key in keys.
    //  Fatal if a zone is not in the mesh. Warns if a zone selects no faces
    selection
    (
        const polyPatch& pp,
        const dictionary& overrides,
        const wordList& keys
    );


    // Member Functions

        //- Zone block names, in dictionary order
        const wordList& names() const
        {
            return names_;
        }

        //- True if zone block j sets keyName
        bool sets(const label j, const word& keyName) const;

        //- True if any zone block sets keyName
        bool setsAny(const word& keyName) const;

        //- Per patch face, the index j of the last zone block that sets
        //  keyName and selects the face, or -1. Size is the patch size
        labelList lastZone(const word& keyName) const;

        //- Warn once per patch if a face is selected by more than one zone
        //  block for the same key
        void overlapWarning() const;
};

}
}
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
