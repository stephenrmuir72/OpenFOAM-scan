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

Class
    Foam::turbulentEpsilonOmegaInletFvPatchScalarField

Group
    grpRASBoundaryConditions grpInletBoundaryConditions

Description
    Inlet boundary condition for the turbulent kinetic energy dissipation
    rate, \c epsilon, or the specific dissipation rate, \c omega. The field
    type is selected from the dimensions of the field to which the condition
    is applied.

    The patch values are calculated using:

    - mixingLength method:

        \f[
            \epsilon_p = \frac{C_{\mu}^{0.75} k^{1.5}}{L}
            \qquad
            \omega_p = \frac{k^{0.5}}{C_{\mu}^{0.25} L}
        \f]

    - viscosityRatio method:

        \f[
            \epsilon_p = \frac{C_{\mu} k^2}{R \nu}
            \qquad
            \omega_p = \frac{k}{R \nu}
        \f]

    where
    \vartable
      \epsilon_p | Patch epsilon values                 [m2/s3]
      \omega_p   | Patch omega values                   [1/s]
      C_\mu      | Empirical model constant
      k          | Turbulent kinetic energy             [m2/s2]
      L          | Mixing length scale                  [m]
      R          | Nominal viscosity ratio, \f$\nu_t/\nu\f$ [-]
      \nu        | Laminar kinematic viscosity on patch [m2/s]
    \endvartable

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries (unmodifiable)
        type            turbulentEpsilonOmegaInlet;

        // Mandatory entries (runtime modifiable)
        method          mixingLength;   // or viscosityRatio
        mixingLength    0.005;          // used when method is mixingLength
        // viscosityRatio 10;           // used when method is viscosityRatio

        // Optional entries (runtime modifiable)
        Cmu             0.09;
        k               k;
        phi             phi;

        // Placeholder
        value           uniform 200;
    }
    \endverbatim

    where the entries mean:
    \table
      Property       | Description                            | Type   | Req'd | Dflt
      method         | Input method: mixingLength or viscosityRatio | word | yes | -
      mixingLength   | Mixing length scale [m]                | scalar | method | -
      viscosityRatio | Nominal viscosity ratio [-]            | scalar | method | -
      Cmu            | Empirical model constant               | scalar |  no   | 0.09
      phi            | Name of flux field                     | word   |  no   | phi
      k              | Name of turbulent kinetic energy field | word   |  no   | k
    \endtable

Note
    - The boundary condition is derived from \c inletOutlet condition.
      Therefore, in the event of reverse flow, a zero-gradient condition
      is applied.
    - The viscosityRatio method assumes the standard eddy-viscosity relations
      \f$\nu_t = C_\mu k^2/\epsilon\f$ and \f$\nu_t = k/\omega\f$.
      The specified ratio is nominal: model limiters or damping, such as the
      SST viscosity limiter, can make the actual \f$\nu_t/\nu\f$ differ from
      the specified value.
    - The order of precedence for the empirical model constant \c Cmu is:
      the \c Cmu entry in the turbulence-model coefficient dictionary,
      the boundary condition dictionary, and the default value of 0.09.
    - The field must have dimensions of either epsilon [m2/s3] or omega [1/s],
      otherwise a fatal error is raised.
    - The empirical model constant \c Cmu is not a spatiotemporal variant field.
      Therefore, the use of the boundary condition may not be fully consistent
      with the turbulence models where \c Cmu is a variant field, such as
      \c realizableKE closure model in this respect.

See also
    Foam::inletOutletFvPatchField

SourceFiles
    turbulentEpsilonOmegaInletFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentEpsilonOmegaInletFvPatchScalarField_H
#define turbulentEpsilonOmegaInletFvPatchScalarField_H

#include "inletOutletFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class turbulentEpsilonOmegaInletFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class turbulentEpsilonOmegaInletFvPatchScalarField
:
    public inletOutletFvPatchScalarField
{
    typedef turbulentEpsilonOmegaInletFvPatchScalarField this_bctype;
    typedef inletOutletFvPatchScalarField parent_bctype;

    // Private Data

        //- Name of the turbulent kinetic energy field
        word kName_;

        //- Use the viscosity ratio method instead of the mixing length method
        bool viscosityRatioMethod_;

        //- True if the field is epsilon, false if it is omega
        bool epsilon_;

        //- Mixing length scale
        scalar mixingLength_;

        //- Viscosity ratio, nu_t/nu
        scalar viscosityRatio_;

        //- Empirical model constant
        scalar Cmu_;


protected:

    // Protected Constructors

        //- Construct from patch, internal field, dictionary and method name
        turbulentEpsilonOmegaInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&,
            const word& method
        );


public:

    //- Runtime type information
    TypeName("turbulentEpsilonOmegaInlet");


    // Constructors

        //- Construct from patch and internal field
        turbulentEpsilonOmegaInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentEpsilonOmegaInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        turbulentEpsilonOmegaInletFvPatchScalarField
        (
            const this_bctype&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        turbulentEpsilonOmegaInletFvPatchScalarField
        (
            const this_bctype&,
            const DimensionedField<scalar, volMesh>&
        );

        //- No copy without an internal field
        turbulentEpsilonOmegaInletFvPatchScalarField
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


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        //- Inherit assignment
        using parent_bctype::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
