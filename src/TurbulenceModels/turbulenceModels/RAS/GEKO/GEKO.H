/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Keysight
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
    \verbatim
        GEKOCoeffs
        {
            Cmu                 0.09;
            C1                  1.44;
            C2                  1.92;
            C3                  -0.33;
            sigmak              1.0;
            sigmaEps            1.3;
            twoLayerTreatment   false;
        }
    \endverbatim

SourceFiles
    GEKO.cxx

\*---------------------------------------------------------------------------*/

#ifndef Foam_GEKO_H
#define Foam_GEKO_H

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

    // Operand Fields

        //- Wall-distance field [m]
        const volScalarField& y_;

        //- Turbulent kinetic energy field [m^2/s^2]
        volScalarField k_;

        //- Specific dissipation rate field [1/s]
        volScalarField omega_;


    // GEKO User-Coefficient Fields

        //- Field to optimise flow separation from smooth surfaces such as
        //- boundary layers etc.
        volScalarField CSEP_;

        //- Field to optimise flows in non-equilibrium near-wall regions
        //- for such as heat transfer, skin friction etc.
        volScalarField CNW_;

        //- Field to optimise spreading rates for free flows.
        volScalarField CMIX_;

        //- Field to optimise round jet flows.
        volScalarField CJET_;


    // GEKO Base Constants (MM:Appendix A, Tables 3-4)

        dimensionedScalar Cmu_;             // 0.09
        dimensionedScalar Cw2_;             // 0.083
        dimensionedScalar sigmaKtilde_;     // 1.0
        dimensionedScalar sigmaOmegaTilde_; // 1.17
        dimensionedScalar kappa_;           // 0.41
        dimensionedScalar CFbTurb_;         // 2.0
        dimensionedScalar CFbLam_;          // 1.0 (25.0 with transition)
        dimensionedScalar A1_;              // 15
        dimensionedScalar Cc1_;             // 1.7
        dimensionedScalar Cc2_;             // 1.4


    // Model-Mode Properties

        //- Flag to enable wall-distance-free variant
        Switch wallDistanceFree_;  // false

        //- Flag to disable Kato-Launder correction for production term in
        //- stagnation regions (MM:p. 5)
        Switch useKatoLaunder_;    // true

        //- Production limiter for turbulent kinetic energy term
        dimensionedScalar CPKlim_;          // 10.0

        //- Realizability limiter constant
        dimensionedScalar CRealize_;        // 0.577


    // Private Member Functions

    // Convenience Functions

        //- Square of the magnitude of the strain-rate tensor field
        tmp<volScalarField> S2() const;

        //- Magnitude of the symmetric strain-rate tensor field
        tmp<volScalarField> Smag() const;

        //- Magnitude of the skew-symmetric rotation-rate tensor field
        tmp<volScalarField> Omag() const;


    // GEKO Helper Functions

        //- Near-wall damping function (MM:Eq. A10)
        tmp<volScalarField> fD() const;

        //- Blending function (MM:Eq. A8)
        tmp<volScalarField> FBlend() const;

        //- Separation calibration function (MM:Eq. A9)
        tmp<volScalarField> FSEP(const volScalarField& Fblend) const;

        //- Jet function (MM:Eq. A14)
        tmp<volScalarField> FJET
        (
            const volScalarField& Sm,
            const volScalarField& Om
        ) const;

        //- Mixing function multiplying the destruction term in omega equation
        //- (MM:Eq. A13)
        tmp<volScalarField> FMIX
        (
            const volScalarField& Fsep,
            const volScalarField& Fblend,
            const volScalarField& Fjet
        ) const;

        //- Near-wall cross-diffusion modifier (MM:Eq. A11)
        tmp<volScalarField> FNW
        (
            const volScalarField& CD,
            const volScalarField& Fd,
            const volScalarField& Destr
        ) const;

        //- Diffusion coefficient for k (MM:Eq. A17)
        tmp<volScalarField> sigmaK(const volScalarField& Fblend) const;

        //- Diffusion coefficient for omega (MM:Eq. A17)
        tmp<volScalarField> sigmaOmega(const volScalarField& Fblend) const;

        //- Coefficient field multiplying Pw in omega equation (MM:Eq. A18)
        tmp<volScalarField> Cw1(const volScalarField& sigmaOmega) const;

        //- Cross-diffusion term field (MM:Eq. A6)
        tmp<volScalarField> CDterm(const volScalarField& sigmaOmega) const;

        //- Production limiter (MM: A15)
        tmp<volScalarField> limitedPk(const volScalarField& Pk) const;

        //- Apply the realizability limiter to nut (MM:Eq. A16)
        void realizabilityLimitNut(volScalarField& nut);


    // Virtual Functions

        //- Calculate the turbulent viscosity
        virtual void correctNut();

        //- The k-source term being used by finite-volume options
        virtual tmp<fvScalarMatrix> kSource() const;

        //- The omega-source term being used by finite-volume options
        virtual tmp<fvScalarMatrix> omegaSource() const;


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

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const { return k_; }

        //- Return the specific dissipation rate
        virtual tmp<volScalarField> omega() const { return omega_; }


    // Evaluation

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>::New
            (
                "DkEff",
                (this->nut_/sigmaKtilde_ + this->nu())  // check the sigma
            );
        }

        //- Return the effective diffusivity for omega
        tmp<volScalarField> DomegaEff() const
        {
            return tmp<volScalarField>::New
            (
                "DomegaEff",
                (this->nut_/sigmaOmegaTilde_ + this->nu())  // check the sigma
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
    #include "GEKO.cxx"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
