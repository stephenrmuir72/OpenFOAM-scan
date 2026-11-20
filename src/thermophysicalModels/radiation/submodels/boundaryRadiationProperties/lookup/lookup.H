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

Class
    Foam::radiation::lookup

Description
    Lookup type of boundary radiation properties.

    The emissivity is a mandatory scalar. The absorptivity is an optional
    scalar; it is an error only when absorptivity is requested for a face
    that has no value. The transmissivity is an optional scalar (default 0).
    Dictionary values such as mappedFile are not accepted. The values are
    read once at construction.

    Optional faceZoneOverrides blocks set emissivity, absorptivity and/or
    transmissivity on the faces of a faceZone. For a face, the last zone in
    dictionary order that sets a property wins. Faces not selected keep the
    patch values.

Usage
    \verbatim

        <patchName>
        {
            type            lookup;
            absorptivity    0.1;
            emissivity      0.1;
            transmissivity  0;

            faceZoneOverrides
            {
                <faceZoneName>
                {
                    emissivity  0.8;
                }
            }
        }

    \endverbatim


\*---------------------------------------------------------------------------*/

#ifndef Foam_radiation_lookup_H
#define Foam_radiation_lookup_H

#include "boundaryRadiationPropertiesPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                            Class lookup Declaration
\*---------------------------------------------------------------------------*/

class lookup
:
    public boundaryRadiationPropertiesPatch
{
protected:

    // Protected data

        //- Reference to the polyPatch
        const polyPatch& pp_;

        //- Emissivity, absorptivity and transmissivity of the patch
        scalar e0_;
        scalar a0_;
        scalar t0_;

        //- True if absorptivity is set in the patch entry
        bool aSet_;

        //- Number of patch faces with no absorptivity value
        label nAMissing_;

        //- Per-face values, only if a faceZone block overrides the patch
        //  value. Empty otherwise, so the scalar is returned for every face
        scalarField eFace_;
        scalarField aFace_;
        scalarField tFace_;

        //- Per-face flag: true if aFace_ holds a value for the face
        boolList aFaceSet_;


public:

    //- Runtime type information
    TypeName("lookup");

    // Constructors

        //- Construct from components
        lookup(const dictionary& dict, const polyPatch& pp);


    //- Destructor
    virtual ~lookup() = default;


    // Member Functions

        //- Return emissivity
        virtual tmp<scalarField> e
        (
            const label bandI,
            const vectorField* incomingDirection,
            const scalarField* T
        ) const;

        //- Return emissivity on face
        virtual scalar e
        (
            const label faceI,
            const label bandI,
            const vector& dir,
            const scalar T
        ) const;

        //- Return absorptivity on patch. Fatal if a face has no value
        virtual tmp<scalarField> a
        (
            const label bandI,
            const vectorField* incomingDirection,
            const scalarField* T
        ) const;

        //- Return absorptivity on face
        virtual scalar a
        (
            const label faceI,
            const label bandI,
            const vector& dir,
            const scalar T
        ) const;

        //- Return transmissivity on patch (default: 0)
        virtual tmp<scalarField> t
        (
            const label bandI,
            const vectorField* incomingDirection,
            const scalarField* T
        ) const;

        //- Return transmissivity on face (default: 0)
        virtual scalar t
        (
            const label faceI,
            const label bandI,
            const vector& dir,
            const scalar T
        ) const;

        //- Return specular reflectivity on patch
        virtual tmp<scalarField> rSpec
        (
            const label bandI,
            const vectorField* incomingDirection,
            const scalarField* T
        ) const;

        //- Return specular reflectivity on face
        virtual scalar rSpec
        (
            const label faceI,
            const label bandI,
            const vector& dir,
            const scalar T
        ) const;

        //- Return diffusive reflectivity on patch
        virtual tmp<scalarField> rDiff
        (
            const label bandI,
            const vectorField* incomingDirection,
            const scalarField* T
        ) const;

        //- Return diffusive reflectivity on face
        virtual scalar rDiff
        (
            const label faceI,
            const label bandI,
            const vector& dir,
            const scalar T
        ) const;

        //- Is Grey
        virtual bool isGrey() const;

        //- Number of bands
        virtual label nBands() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
