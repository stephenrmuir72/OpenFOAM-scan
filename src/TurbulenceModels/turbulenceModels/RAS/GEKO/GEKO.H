/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Keysight Ltd.
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
    Foam::RASModels::GEKO

Group
    grpRASTurbulence

Description


    Reference:
    \verbatim
        Standard model (tag:MM):
            Menter, F. R., & Matyushenko, A. (2025).
            Generalized k−omega (GEKO) Two-Equation Turbulence Model.
            AIAA Journal, 63(11), 4590-4606.
    \endverbatim

Usage

    The default model coefficients are

SourceFiles
    GEKO.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_RASModels_GEKO_H
#define Foam_RASModels_GEKO_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class GEKO Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class GEKO
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{
    // Private Data

    // Model switches

        //- Flag to enable wall-distance-free model variant
        Switch wallDistanceFree_;      // (MM:p. 14; if true: F_blend = 1)

        //- Flag to disable production limiter for turbulent kinetic energy term
        Switch productionLimiter_;     // (MM:Eq. 16; default true)

        //- Flag to disable Kato-Launder correction for production term in
        //- stagnation regions
        Switch katoLaunder_;           // (MM:Eq. 15)

        //- Flag to enable dilatation correction term in the Reynolds stress
        //- tensor
        Switch dilatationCorrection_;  // (MM:Eq. A7)


    // Model coefficients

        // (MM:Table 3; Main model coefficients for GEKO model)
        dimensionedScalar Cmu_;
        dimensionedScalar Comega2_;
        dimensionedScalar sigmakTilde_;
        dimensionedScalar sigmaOmegaTilde_;
        dimensionedScalar CRealize_;

        // (MM:Table 4; Auxiliary model coefficients for GEKO model)
        dimensionedScalar kappa_;
        dimensionedScalar CFbTurb_;
        dimensionedScalar CFbLam_;
        dimensionedScalar Aplus_;
        dimensionedScalar Cc1_;
        dimensionedScalar Cc2_;

        // (MM:Eq. 16)
        dimensionedScalar CPKlim_;


    // Operand Fields

        //- Wall-distance field [m]
        const volScalarField& y_;

        //- Turbulent kinetic energy field [m^2/s^2]
        volScalarField k_;

        //- Specific dissipation rate field [1/s]
        volScalarField omega_;


    // GEKO model fields

        //- Field to optimise flow separation from smooth surfaces such as
        //- boundary layers etc.
        volScalarField CSEP_;

        //- Field to optimise flows in non-equilibrium near-wall regions
        //- for such as heat transfer, skin friction etc.
        volScalarField CNW_;

        //- Field to optimise spreading rates for free flows
        volScalarField CMIX_;

        //- Field to optimise round jet flows
        volScalarField CJET_;


    // Auxillary fields

        //- Blending function field (MM:Eq. A8)
        volScalarField FBlend_;


    // Optional machine-learning fields (MM:p. 13)

        //- Machine-learning field to optimise turbulent kinetic energy
        volScalarField Ck_;

        //- Machine-learning field to optimise specific dissipation rate
        volScalarField Comega_;


    // Private Member Functions

    // Convenience Functions



    // GEKO Helper Functions


    // Virtual Functions

        //- Calculate the turbulent viscosity
        virtual void correctNut();

public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("GEKO");


    // Generated Methods

        //- No copy construct
        GEKO(const GEKO&) = delete;

        //- No copy assignment
        void operator=(const GEKO&) = delete;


    // Constructors

        //- Construct from components
        GEKO
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~GEKO() = default;


    // Member Functions

    // Access

        //- Return the turbulence kinetic energy field [m^2/s^2]
        virtual tmp<volScalarField> k() const { return k_; }

        //- Return the specific dissipation rate field [1/s]
        virtual tmp<volScalarField> omega() const { return omega_; }


    // Evaluation

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>::New
            (
                "DkEff",
                (this->nut_/sigmaKtilde_ + this->nu())
            );
        }

        //- Return the effective diffusivity for omega
        tmp<volScalarField> DomegaEff() const
        {
            return tmp<volScalarField>::New
            (
                "DomegaEff",
                (this->nut_/sigmaOmegaTilde_ + this->nu())
            );
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();


    // I-O

        //- Re-read model coefficients if they have changed
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "GEKO.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
