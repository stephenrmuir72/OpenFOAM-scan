/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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
    Foam::ignoreFaPatchField

Description
    A boundary condition for invalid edges which should be ignored.
    Currently functionally identical to a zero-gradient condition.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type    ignore;
    }
    \endverbatim

SourceFiles
    ignoreFaPatchField.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_ignoreFaPatchField_H
#define Foam_ignoreFaPatchField_H

#include "zeroGradientFvPatchField.H"
#include "areaFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class ignoreFaPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class ignoreFaPatchField
:
    public zeroGradientFvPatchField<Type>
{
    typedef ignoreFaPatchField<Type> this_bctype;
    typedef zeroGradientFvPatchField<Type> parent_bctype;

public:

    //- Runtime type information
    TypeName("ignore");


    // Constructors

        //- Construct from patch and internal field
        ignoreFaPatchField
        (
            const faPatch& p,
            const DimensionedField<Type, areaMesh>& iF
        )
        :
            parent_bctype(p, iF)
        {}

        //- Construct from patch, internal field and dictionary
        ignoreFaPatchField
        (
            const faPatch& p,
            const DimensionedField<Type, areaMesh>& iF,
            const dictionary& dict
        )
        :
            parent_bctype(p, iF, dict)
        {}

        //- Construct by mapping onto a new patch
        ignoreFaPatchField
        (
            const ignoreFaPatchField<Type>& ptf,
            const faPatch& p,
            const DimensionedField<Type, areaMesh>& iF,
            const faPatchFieldMapper& mapper
        )
        :
            parent_bctype(ptf, p, iF, mapper)
        {}

        //- Copy construct
        ignoreFaPatchField(const this_bctype& ptf)
        #ifdef Foam_fvPatchField_copyConstruct
        : ignoreFaPatchField<Type>(ptf, ptf.internalField()) {}
        #else
        = delete;
        #endif

        //- Copy construct with internal field reference
        ignoreFaPatchField
        (
            const ignoreFaPatchField<Type>& ptf,
            const DimensionedField<Type, areaMesh>& iF
        )
        :
            parent_bctype(ptf, iF)
        {}

        //- Clone with an internal field reference
        virtual tmp<faPatchField<Type>> clone
        (
            const DimensionedField<Type, areaMesh>& iF
        ) const
        {
            return faPatchField<Type>::Clone(*this, iF);
        }


    //- Destructor
    virtual ~ignoreFaPatchField() = default;


    // Member Operators

        //- Inherit assignment
        using parent_bctype::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
