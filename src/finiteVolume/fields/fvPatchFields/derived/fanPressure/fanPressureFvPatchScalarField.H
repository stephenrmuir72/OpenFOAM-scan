/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd
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
    Foam::fanPressureFvPatchScalarField

Group
    grpInletBoundaryConditions grpOutletBoundaryConditions

Description
    This boundary condition can be applied to assign either a pressure inlet
    or outlet total pressure condition for a fan.

    The pressure drop \c deltaP can be optionally under-relaxed using the
    \c relax factor (0 < relax <= 1). The update is:

        deltaP = relax*deltaPNew + (1 - relax)*deltaP0

    where \c deltaP0 is the reference value from the previous update. The
    reference is refreshed only when the time index changes, so the
    relaxation acts within the iterations of a time step and is therefore
    also valid for transient simulations.

    When \c relax is active, \c deltaP0 is written and read automatically to
    preserve the relaxation state across restarts; it is not intended as a
    normal user configuration entry.

Usage
    \table
        Property     | Description                       | Required | Default
        fanCurve     | Pressure vs flow-rate             | yes |
        direction    | direction of flow through fan [in/out] | yes |
        p0           | environmental total pressure      | yes |
        nonDimensional | uses non-dimensional table      | no  | false
        relax         | under-relaxation factor for deltaP    | no | -1
        file         | fan curve file name               | legacy |
        outOfBounds  | out of bounds handling            | legacy |
    \endtable

    For the inputs and formulation related to non-dimensional operation, refer
    to \link fanNonDimensional.H \endlink.

    The inherited entries are elaborated in:
      - \link totalPressureFvPatchScalarField.H \endlink

    Example of the boundary condition specification:
    \verbatim
    inlet
    {
        type            fanPressure;
        direction       in;
        fanCurve
        {
            type table;
            file "<constant>/fanCurve";
            outOfBounds clamp;  // Optional out-of-bounds handling
        }
        p0              uniform 0;
        value           uniform 0;
    }

    // Legacy specification
    outlet
    {
        type            fanPressure;
        direction       out;
        file            "<constant>/fanCurve";
        outOfBounds     clamp;
        p0              uniform 0;
        value           uniform 0;
    }
    \endverbatim

Note
    For compatibility with older versions (OpenFOAM-v2006 and earlier),
    a missing \c fanCurve keyword is treated as a tableFile and makes
    the \c file keyword mandatory.

See also
    Foam::fanNonDimensional
    Foam::fanFvPatchField
    Foam::totalPressureFvPatchScalarField
    Foam::Function1

SourceFiles
    fanPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fanPressureFvPatchScalarField_H
#define Foam_fanPressureFvPatchScalarField_H

#include "totalPressureFvPatchScalarField.H"
#include "fanNonDimensional.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class fanPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class fanPressureFvPatchScalarField
:
    public totalPressureFvPatchScalarField
{
    typedef fanPressureFvPatchScalarField this_bctype;
    typedef totalPressureFvPatchScalarField parent_bctype;

public:

        //- Fan flow direction
        enum fanFlowDirection : uint8_t
        {
            ffdIn,
            ffdOut
        };

        //- Fan flow direction names
        static const Enum<fanFlowDirection> fanFlowDirectionNames_;


private:

    // Private Data

        //- Run-time selectable fan curve
        autoPtr<Function1<scalar>> fanCurve_;

        //- Use non-dimensional curve
        bool nonDimensional_;

        //- Non-dimensional fan-curve formulation
        fanNonDimensional fanNonDimensional_;

        //- Under-relaxation factor for the pressure drop
        scalar relaxFactor_;

        //- Reference pressure drop for under-relaxation
        scalar pdFan0_;

        //- Time index of the reference pressure drop
        label timeIndex_;

        //- Direction of flow through the fan relative to patch
        fanFlowDirection direction_;

public:

    //- Runtime type information
    TypeName("fanPressure");


    // Constructors

        //- Construct from patch and internal field
        fanPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fanPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        fanPressureFvPatchScalarField
        (
            const this_bctype&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        fanPressureFvPatchScalarField
        (
            const this_bctype&,
            const DimensionedField<scalar, volMesh>&
        );

        //- No copy without an internal field
        fanPressureFvPatchScalarField(const this_bctype&) = delete;

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return fvPatchField<scalar>::Clone(*this, iF);
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream& os) const;


    // Member Operators

        //- Inherit assignment
        using parent_bctype::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
