/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 OpenCFD Ltd.
    Copyright (C) Keysight Technologies
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
    Foam::meshPointPatchDisplacementVectorField

Description
    Constrain the field according to the point constraints (on meshPointPatch).
    Can only be used on meshPointPatch and displacementMotionSolver
    (to provide point0).

    Action depends on pointConstraint (stored on meshPointPatch) and the
    fixedValue (from the pointPatchField).
    - 2 (surface) : attract to nearest point on plane. Plane indicated
                    by normal (from pointConstraint) and base point
                    (given as points0 + fixedValue).
    - 1 (edge)    : attract to nearest point on edge. Edge indicated
                    by direction (from pointConstraint) and base point
                    (given as points0 + fixedValue).
    - 0 (point)   : frozen to base point.

SourceFiles
    meshPointPatchDisplacementVectorField.cxx

\*---------------------------------------------------------------------------*/

#ifndef meshPointPatchDisplacementVectorField_H
#define meshPointPatchDisplacementVectorField_H

#include "fixedValuePointPatchField.H"
#include "meshPointPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class meshPointPatchDisplacementVectorField Declaration
\*---------------------------------------------------------------------------*/

class meshPointPatchDisplacementVectorField
:
    public fixedValuePointPatchField<vector>
{
private:

    // Private data

        //- Local reference cast into the meshPoint patch
        const meshPointPatch& mpPatch_;


        //- No copy assignment
        void operator=
        (
            const meshPointPatchDisplacementVectorField&
        ) = delete;


public:

    //- Runtime type information
    TypeName("displacementMeshPoint");


    // Constructors

        //- Construct from patch and internal field
        meshPointPatchDisplacementVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        meshPointPatchDisplacementVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patch field onto a new patch
        meshPointPatchDisplacementVectorField
        (
            const meshPointPatchDisplacementVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<vector>> clone() const
        {
            return autoPtr<pointPatchField<vector>>
            (
                new meshPointPatchDisplacementVectorField
                (
                    *this,
                    internalField()
                )
            );
        }

        //- Construct as copy setting internal field reference
        meshPointPatchDisplacementVectorField
        (
            const meshPointPatchDisplacementVectorField&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<vector>> clone
        (
            const DimensionedField<vector, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<vector>>
            (
                new meshPointPatchDisplacementVectorField
                (
                    *this,
                    iF
                )
            );
        }

    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Evaluate the patch field, sets updated() to false
        virtual void evaluate
        (
            const Pstream::commsTypes commsType = Pstream::commsTypes::buffered
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
