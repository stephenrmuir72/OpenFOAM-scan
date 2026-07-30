/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2026 OpenCFD Ltd.
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
    Foam::edgeNormalFixedValueFaPatchVectorField

Description
    Edge-normal fixed value vector finite-area boundary condition.

    Describes a surface-normal vector boundary condition by its magnitude.
    Note: The value is positive for outward-pointing vectors.

    In the finiteArea refactor this is a thin backward-compatibility alias
    that preserves the legacy dictionary type-name "edgeNormalFixedValue".
    All numerics are inherited unchanged from the finite-volume
    \c surfaceNormalFixedValue boundary condition: both store a scalar
    reference magnitude and set the patch value to
    \c refValue*(outward patch normal). On the reduced-order fvMesh the
    finite-area edge normal equals the finite-volume face normal
    (\c patch().nf() == Sf/|Sf| == the legacy \c edgeNormals()), so the two
    are numerically identical. New cases should prefer
    \c surfaceNormalFixedValue directly.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            edgeNormalFixedValue;
        refValue        uniform -1;   // 1 INTO the domain
    }
    \endverbatim

See also
    Foam::surfaceNormalFixedValueFvPatchVectorField

SourceFiles
    edgeNormalFixedValueFaPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_edgeNormalFixedValueFaPatchVectorField_H
#define Foam_edgeNormalFixedValueFaPatchVectorField_H

#include "surfaceNormalFixedValueFvPatchVectorField.H"
#include "areaFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class edgeNormalFixedValueFaPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class edgeNormalFixedValueFaPatchVectorField
:
    public surfaceNormalFixedValueFvPatchVectorField
{
    typedef edgeNormalFixedValueFaPatchVectorField this_bctype;
    typedef surfaceNormalFixedValueFvPatchVectorField parent_bctype;

public:

    //- Runtime type information
    TypeName("edgeNormalFixedValue");


    // Constructors

        //- Construct from patch and internal field
        edgeNormalFixedValueFaPatchVectorField
        (
            const faPatch&,
            const DimensionedField<vector, areaMesh>&
        );

        //- Construct from patch, internal field and dictionary
        edgeNormalFixedValueFaPatchVectorField
        (
            const faPatch&,
            const DimensionedField<vector, areaMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        edgeNormalFixedValueFaPatchVectorField
        (
            const this_bctype&,
            const faPatch&,
            const DimensionedField<vector, areaMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        edgeNormalFixedValueFaPatchVectorField
        (
            const this_bctype&,
            const DimensionedField<vector, areaMesh>&
        );

        //- Copy construct
        edgeNormalFixedValueFaPatchVectorField(const this_bctype& ptf)
        #ifdef Foam_fvPatchField_copyConstruct
        : this_bctype(ptf, ptf.internalField()) {}
        #else
        = delete;
        #endif

        //- Clone with an internal field reference
        virtual tmp<faPatchField<vector>> clone
        (
            const DimensionedField<vector, areaMesh>& iF
        ) const
        {
            return faPatchField<vector>::Clone(*this, iF);
        }


    // Member Operators

        //- Inherit assignment
        using parent_bctype::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
