/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    GEKO (Generalized k-omega) turbulence model base class.

    Implements the GEKO model as defined in:
      Menter, F.R. & Matyushenko, A. (2025)
      "Generalized k-ω turbulence model for RANS simulations ..."
    using Appendix A equations (A1)-(A18) and related notes.

    - Compressible form (ρ-weighted transport)
    - Production limiter (Eq. 16 / A15) and realizability limiter (Eq. 17 / A16)
    - Optional Kato–Launder production (Eq. 15)
    - Optional compressible dilatation correction in τ_ij (note near A7)
    - Coefficients CSEP, CNW, CMIX, CJET are implemented as volScalarField

\*---------------------------------------------------------------------------*/

#ifndef GEKOBase_H
#define GEKOBase_H

#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class BasicEddyViscosityModel>
class GEKOBase
:
    public BasicEddyViscosityModel
{
protected:

    // Typedefs from the parent model (same pattern as kOmegaSSTBase)
    typedef typename BasicEddyViscosityModel::alphaField alphaField;
    typedef typename BasicEddyViscosityModel::rhoField rhoField;
    typedef typename BasicEddyViscosityModel::transportModel transportModel;

    // ---- Model constants (Appendix A Tables 3/4, Eq. 16, Eq. 17) ----

    // Table 3
    dimensionedScalar Cmu_;             // C_mu
    dimensionedScalar Comega2_;         // C_omega2
    dimensionedScalar sigmakTilde_;     // \tilde\sigma_k
    dimensionedScalar sigmaomegaTilde_; // \tilde\sigma_omega
    dimensionedScalar CRealize_;        // C_Realize (= 1/sqrt(3))

    // Table 4 + Appendix A defs
    dimensionedScalar kappa_;           // von Karman constant
    dimensionedScalar CFbTurb_;         // C_FbTurb
    dimensionedScalar CFbLam_;          // C_FbLam
    dimensionedScalar APlus_;           // A^+
    dimensionedScalar Cc1_;             // C_c1
    dimensionedScalar Cc2_;             // C_c2

    // Eq. 16
    dimensionedScalar CPKlim_;          // C_PKlim

    // ---- Switches (per your A1 answers) ----
    Switch wallDistanceFree_;           // if yes: F_blend = 1 everywhere
    Switch productionLimiter_;          // default yes
    Switch katoLaunder_;                // optional Eq. 15
    Switch dilatationCorrection_;       // option for τ_ij compressible term (note near A7)

    // ---- Fields ----
    const volScalarField& y_;           // wall distance (from wallDist)

    volScalarField k_;
    volScalarField omega_;

    // GEKO coefficient fields (your requirement A1.5)
    volScalarField CSEP_;
    volScalarField CNW_;
    volScalarField CMIX_;
    volScalarField CJET_;

    // Optional ML coefficient fields (Appendix note near A7)
    volScalarField Ck_;
    volScalarField Comega_;


    // ---- Internal helpers ----
    virtual void correctNut(const volScalarField& S);
    virtual void correctNut();

public:

    //- Runtime type information is provided by the wrapper class (RASModels::GEKO)

    // Constructors
    GEKOBase
    (
        const word& type,
        const alphaField& alpha,
        const rhoField& rho,
        const volVectorField& U,
        const surfaceScalarField& alphaRhoPhi,
        const surfaceScalarField& phi,
        const transportModel& transport,
        const word& propertiesName
    );

    //- Destructor
    virtual ~GEKOBase() = default;


    // Member Functions

    //- Read model coefficients (from GEKOCoeffs)
    virtual bool read();

    //- Return turbulence kinetic energy [m^2/s^2]
    virtual tmp<volScalarField> k() const
    {
        return k_;
    }

    //- Return specific dissipation rate [1/s]
    virtual tmp<volScalarField> omega() const
    {
        return omega_;
    }

    //- Return turbulence dissipation rate epsilon = Cmu*k*omega
    virtual tmp<volScalarField> epsilon() const;

    //- Correct the turbulence model
    virtual void correct();

    //- Reynolds stress tensor (kinematic) for output; includes optional dilatation correction term
    virtual tmp<volSymmTensorField> R() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
