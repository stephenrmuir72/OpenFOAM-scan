/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::turbulentMixingLengthFrequencyInletFvPatchScalarField

Group
    grpRASBoundaryConditions grpInletBoundaryConditions

Description
    This boundary condition provides a turbulence specific dissipation,
    \f$\omega\f$ (omega) inlet condition based on a specified mixing length.
    The patch values are calculated using:

        \f[
            \omega_p = \frac{k^{0.5}}{C_{\mu}^{0.25} L}
        \f]

    where

    \vartable
        \omega_p | patch omega values
        C_{\mu} | Empirical model constant (default 0.09)
        k       | turbulence kinetic energy
        L       | length scale
    \endvartable

Usage
    \table
      Property     | Description                  | Type   | Req'd | Dflt
      mixingLength | Mixing length scale [m]      | scalar |  yes  | -
      Cmu          | Empirical model constant     | scalar |  no   | 0.09
      phi          | Name of flux field           | word   |  no   | phi
      k            | Name of turbulent kinetic energy field | word | no | k
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            turbulentMixingLengthFrequencyInlet;
        mixingLength    0.005;
        Cmu             0.09;          // optional
        value           uniform 200;   // placeholder
    }
    \endverbatim

Note
    - In the event of reverse flow, a zero-gradient condition is applied.
    - The order of precedence for the empirical model constant \c Cmu is:
      the \c Cmu entry in the turbulence-model coefficient dictionary,
      the boundary condition dictionary, and the default value of 0.09.

See also
    Foam::turbulentEpsilonOmegaInletFvPatchScalarField


SourceFiles
    turbulentMixingLengthFrequencyInletFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentMixingLengthFrequencyInletFvPatchScalarField_H
#define turbulentMixingLengthFrequencyInletFvPatchScalarField_H

#include "turbulentEpsilonOmegaInletFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
    Class turbulentMixingLengthFrequencyInletFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class turbulentMixingLengthFrequencyInletFvPatchScalarField
:
    public turbulentEpsilonOmegaInletFvPatchScalarField
{
    typedef turbulentMixingLengthFrequencyInletFvPatchScalarField this_bctype;
    typedef turbulentEpsilonOmegaInletFvPatchScalarField parent_bctype;


public:

    //- Runtime type information
    TypeName("turbulentMixingLengthFrequencyInlet");


    // Constructors

        //- Construct from patch and internal field
        turbulentMixingLengthFrequencyInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentMixingLengthFrequencyInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        turbulentMixingLengthFrequencyInletFvPatchScalarField
        (
            const this_bctype&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        turbulentMixingLengthFrequencyInletFvPatchScalarField
        (
            const this_bctype&,
            const DimensionedField<scalar, volMesh>&
        );

        //- No copy without an internal field
        turbulentMixingLengthFrequencyInletFvPatchScalarField
        (
            const this_bctype&
        ) = delete;

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return fvPatchField<scalar>::Clone(*this, iF);
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
