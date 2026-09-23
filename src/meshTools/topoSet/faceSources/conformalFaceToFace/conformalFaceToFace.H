/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Patrick Wilhelm, Data Lab Hell
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
    Foam::conformalFaceToFace

Description
    A \c topoSetFaceSource to select the boundary faces of this region that are
    conformal with a faceSet of boundary faces in another region.

    The target is searched for one \c distance along the source face normal,
    and a conformal target opposes that normal. Use \c flipDistance when this
    region is stacked behind the source, or \c anySign to accept either
    stacking. \c tolerance bounds the centre match and \c areaTolerance the
    area difference from the source face.

    Operands:
    \table
      Operand | Type    | Location
      input   | faceSet | $FOAM_CASE/\<region\>/constant/polyMesh/sets/\<set\>
      output  | faceSet | $FOAM_CASE/constant/polyMesh/sets/\<set\>
    \endtable

Usage
    Minimal example by using \c system/topoSetDict.actions:
    \verbatim
    {
        // Mandatory (inherited) entries
        name        <name>;
        type        faceSet;
        action      <action>;

        // Mandatory entries
        source      conformalFaceToFace;
        region      <regionName>;
        set         <faceSetName>;

        // Optional entries
        patches     (<patchName1> <patchName2> ... <patchNameN>);
        patch       <patchName>;
        distance    <scalar>;
        flipDistance <bool>;
        anySign     <bool>;
        tolerance   <scalar>;
        areaTolerance <scalar>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                         | Type | Req'd | Dflt
      name       | Name of faceSet                     | word |  yes  | -
      type       | Type name: faceSet                  | word |  yes  | -
      action     | Action applied on faces - see below | word |  yes  | -
      source     | Source name: conformalFaceToFace    | word |  yes  | -
      region     | Region holding the input faceSet    | word |  yes  | -
      set        | Name of the input faceSet           | word |  yes  | -
      patches    | Restrict matching to target patches | wordRes | no | all
      patch      | Restrict matching to a target patch | wordRe | no | all
      distance   | Separation of matched face centres  | scalar | no | 0
      flipDistance | Seek the target against the normal | bool | no | false
      anySign    | Seek either way, no normal test     | bool | no | false
      tolerance  | Match tolerance as a length         | scalar | no | 0
      areaTolerance | Relative bound on the area diff  | scalar | no | 1e-3
    \endtable

    Options for the \c action entry:
    \verbatim
      new      | Create a new faceSet from selected faces
      add      | Add selected faces into this faceSet
      subtract | Remove selected faces from this faceSet
    \endverbatim

Note
    \c patches takes precedence over \c patch.

See also
    - Foam::topoSetSource
    - Foam::topoSetFaceSource
    - Foam::faceToFace

SourceFiles
    conformalFaceToFace.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_conformalFaceToFace_H
#define Foam_conformalFaceToFace_H

#include "topoSetFaceSource.H"
#include "polyMesh.H"
#include "wordRes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class conformalFaceToFace Declaration
\*---------------------------------------------------------------------------*/

class conformalFaceToFace
:
    public topoSetFaceSource
{
    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

        //- Name of the region holding the input faceSet
        const word sourceRegionName_;

        //- Name of the input faceSet in the source region
        const word sourceSetName_;

        //- Target patches to match against, empty for all
        wordRes selectedPatches_;

        //- Separation between matched face centres
        const scalar distance_;

        //- Match tolerance as a length, zero for the per-face bound
        const scalar tolerance_;

        //- Relative bound on the area difference of a matched pair
        const scalar areaTolerance_;

        //- Whether the target normal must oppose the source normal
        const bool directional_;

        //- Whether the target is sought against the source normal
        const bool flip_;

        //- The source region on this rank, read unregistered
        autoPtr<polyMesh> meshSourcePtr_;


    // Private Member Functions

        //- The centre match bound for a source face of the given area
        scalar faceTolerance(const scalar area) const;

        //- Select the conformal target faces
        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("conformalFaceToFace");


    // Constructors

        //- Construct from dictionary
        conformalFaceToFace(const polyMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~conformalFaceToFace() = default;


    // Member Functions

        virtual void applyToSet
        (
            const topoSetSource::setAction action,
            topoSet& set
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
