/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
\*---------------------------------------------------------------------------*/

#include "GEKOBase.H"

#include "wallDist.H"
#include "fvOptions.H"
#include "bound.H"
#include "calculatedFvPatchField.H"

#include <cmath>

namespace Foam
{

template<class BasicEddyViscosityModel>
GEKOBase<BasicEddyViscosityModel>::GEKOBase
(
    const word& type,
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
:
    BasicEddyViscosityModel
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    // --- Defaults from Appendix A (Tables 3/4 and Eq.16/Eq.17) ---
    Cmu_(dimensioned<scalar>::getOrAddToDict("Cmu", this->coeffDict_, 0.09)),
    Comega2_(dimensioned<scalar>::getOrAddToDict("Comega2", this->coeffDict_, 0.083)),
    sigmakTilde_(dimensioned<scalar>::getOrAddToDict("sigmakTilde", this->coeffDict_, 1.0)),
    sigmaomegaTilde_(dimensioned<scalar>::getOrAddToDict("sigmaomegaTilde", this->coeffDict_, 1.17)),
    CRealize_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "CRealize",
            this->coeffDict_,
            scalar(1.0)/std::sqrt(scalar(3.0))  // Eq. 17 / A16: 1/sqrt(3)
        )
    ),

    kappa_(dimensioned<scalar>::getOrAddToDict("kappa", this->coeffDict_, 0.41)),
    CFbTurb_(dimensioned<scalar>::getOrAddToDict("CFbTurb", this->coeffDict_, 2.0)),
    CFbLam_(dimensioned<scalar>::getOrAddToDict("CFbLam", this->coeffDict_, 1.0)),
    APlus_(dimensioned<scalar>::getOrAddToDict("APlus", this->coeffDict_, 15.0)),
    Cc1_(dimensioned<scalar>::getOrAddToDict("Cc1", this->coeffDict_, 1.7)),
    Cc2_(dimensioned<scalar>::getOrAddToDict("Cc2", this->coeffDict_, 1.4)),

    CPKlim_(dimensioned<scalar>::getOrAddToDict("CPKlim", this->coeffDict_, 10.0)),

    // switches
    wallDistanceFree_(this->coeffDict_.getOrDefault<Switch>("wallDistanceFree", false)),
    productionLimiter_(this->coeffDict_.getOrDefault<Switch>("productionLimiter", true)),
    katoLaunder_(this->coeffDict_.getOrDefault<Switch>("katoLaunder", false)),
    dilatationCorrection_(this->coeffDict_.getOrDefault<Switch>("dilatationCorrection", false)),

    // wall distance
    y_(wallDist::New(this->mesh_).y()),

    // turbulence fields
    k_
    (
        IOobject
        (
            IOobject::groupName("k", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    // GEKO coefficient fields (must be present; your requirement)
    CSEP_
    (
        IOobject
        (
            IOobject::groupName("CSEP", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    CNW_
    (
        IOobject
        (
            IOobject::groupName("CNW", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    CMIX_
    (
        IOobject
        (
            IOobject::groupName("CMIX", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    CJET_
    (
        IOobject
        (
            IOobject::groupName("CJET", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    // ML coefficient fields (Appendix note; default should be uniform 0 in case setup)
    Ck_
    (
        IOobject
        (
            IOobject::groupName("Ck", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    Comega_
    (
        IOobject
        (
            IOobject::groupName("Comega", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);
}


template<class BasicEddyViscosityModel>
bool GEKOBase<BasicEddyViscosityModel>::read()
{
    if (BasicEddyViscosityModel::read())
    {
        // read constants if present
        Cmu_.readIfPresent(this->coeffDict());
        Comega2_.readIfPresent(this->coeffDict());
        sigmakTilde_.readIfPresent(this->coeffDict());
        sigmaomegaTilde_.readIfPresent(this->coeffDict());
        CRealize_.readIfPresent(this->coeffDict());

        kappa_.readIfPresent(this->coeffDict());
        CFbTurb_.readIfPresent(this->coeffDict());
        CFbLam_.readIfPresent(this->coeffDict());
        APlus_.readIfPresent(this->coeffDict());
        Cc1_.readIfPresent(this->coeffDict());
        Cc2_.readIfPresent(this->coeffDict());

        CPKlim_.readIfPresent(this->coeffDict());

        // switches
        wallDistanceFree_ = this->coeffDict().getOrDefault<Switch>("wallDistanceFree", wallDistanceFree_);
        productionLimiter_ = this->coeffDict().getOrDefault<Switch>("productionLimiter", productionLimiter_);
        katoLaunder_ = this->coeffDict().getOrDefault<Switch>("katoLaunder", katoLaunder_);
        dilatationCorrection_ = this->coeffDict().getOrDefault<Switch>("dilatationCorrection", dilatationCorrection_);

        return true;
    }

    return false;
}


template<class BasicEddyViscosityModel>
tmp<volScalarField> GEKOBase<BasicEddyViscosityModel>::epsilon() const
{
    // epsilon = Cmu * k * omega
    return volScalarField::New
    (
        IOobject::groupName("epsilon", this->alphaRhoPhi_.group()),
        IOobject::NO_REGISTER,
        Cmu_*k_*omega_
    );
}


template<class BasicEddyViscosityModel>
void GEKOBase<BasicEddyViscosityModel>::correctNut(const volScalarField& S)
{
    // Appendix A, Eq. (A16) / Eq. (17): nu_t = k / max(omega, S/CRealize)
    this->nut_ = k_/max(omega_, S/CRealize_);

    this->nut_.correctBoundaryConditions();

    fv::options::New(this->mesh_).correct(this->nut_);
}


template<class BasicEddyViscosityModel>
void GEKOBase<BasicEddyViscosityModel>::correctNut()
{
    // S = sqrt(2*Sij*Sij), with Sij = symm(gradU)
    const volScalarField S
    (
        sqrt(2.0*magSqr(symm(fvc::grad(this->U_))))
    );

    correctNut(S);
}


template<class BasicEddyViscosityModel>
tmp<volSymmTensorField> GEKOBase<BasicEddyViscosityModel>::R() const
{
    // Output Reynolds stress (kinematic): R = (2/3)k I - nu_t * (2S_ij - (2/3)divU δ_ij)
    // The dilatation correction (-(2/3)divU term) is optional per your requirement.

    tmp<volScalarField> tk(this->k());

    // Patch type list (same pattern as eddyViscosity::R())
    wordList patchFieldTypes
    (
        this->nut_.boundaryField().size(),
        calculatedFvPatchField<symmTensor>::typeName
    );

    // gradU in OpenFOAM is Hessian form (∂U_j/∂x_i), but symm(...) is invariant.
    tmp<volTensorField> tgradU = fvc::grad(this->U_);

    tmp<volSymmTensorField> tR
    (
        volSymmTensorField::New
        (
            IOobject::groupName("R", this->alphaRhoPhi_.group()),
            IOobject::NO_REGISTER,
            ((2.0/3.0)*I)*tk()
          - this->nut_
           *(
                dilatationCorrection_
              ? devTwoSymm(tgradU())
              : (2.0*symm(tgradU()))
            ),
            patchFieldTypes
        )
    );

    return tR;
}


template<class BasicEddyViscosityModel>
void GEKOBase<BasicEddyViscosityModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;

    fv::options& fvOptions(fv::options::New(this->mesh_));

    // Base class bookkeeping (same pattern as kOmegaSSTBase::correct())
    BasicEddyViscosityModel::correct();

    // gradU, strain-rate magnitude S, and vorticity magnitude Omega (paper uses Ω)
    tmp<volTensorField> tgradU = fvc::grad(U);

    const volScalarField S2(2.0*magSqr(symm(tgradU())));
    const volScalarField S(sqrt(S2));

    // Ω = |curl(U)| = sqrt(2*Ω_ijΩ_ij), with Ω_ij = skew(gradU)
    const volScalarField OmegaMag(sqrt(2.0*magSqr(skew(tgradU()))));

    // --- Reuse SST wall-function handling for omega (exact same sequence) ---
    omega_.boundaryFieldRef().updateCoeffs();
    omega_.boundaryFieldRef().template evaluateCoupled<coupledFvPatch>();


    // ---------------------------------------------------------------------
    // Appendix A auxiliary functions and coefficients
    // ---------------------------------------------------------------------

    // k~ = max(k, C_FbLam * nu * omega)   (A9)
    const volScalarField kTilde(max(k_, CFbLam_*this->nu()*omega_));

    // L_T = sqrt(k~)/(Cmu*omega)         (A9)
    const volScalarField LT(sqrt(kTilde)/(Cmu_*omega_));

    // x_blend = C_FbTurb * L_T / y       (A8)
    const volScalarField xBlend(CFbTurb_*LT/y_);

    // F_blend = tanh(x_blend^4)          (A8)
    volScalarField FBlend(tanh(pow4(xBlend)));

    if (wallDistanceFree_)
    {
        FBlend = scalar(1);
    }

    // sigma_k and sigma_omega             (A17)
    const volScalarField sigmaK
    (
        sigmakTilde_*(scalar(1) + (CSEP_ - scalar(1))*scalar(0.25)*FBlend)
    );

    const volScalarField sigmaOmega
    (
        sigmaomegaTilde_*(scalar(1) + (CSEP_ - scalar(1))*FBlend)
    );

    // C_omega1                              (A18)
    const scalar sqrtCmu = std::sqrt(Cmu_.value());
    const volScalarField Comega1
    (
        (scalar(1)/Cmu_)
       *(
            Comega2_
          - sqr(kappa_)*sqrtCmu/sigmaOmega
        )
    );

    // y~+                                    (A10)
    const volScalarField yTildePlus
    (
        (scalar(1)/kappa_)*(k_/omega_)/this->nu()
    );

    // f_D                                    (A10)
    const volScalarField fD
    (
        scalar(1)/(scalar(1) + sqr(yTildePlus/APlus_))
    );

    // C_psi = 0.2454 * CSEP^{-0.803}         (A9)
    // Implemented via explicit loops to avoid relying on field-pow overloads.
    volScalarField Cpsi(scalar(0)*FBlend);

    forAll(Cpsi.primitiveField(), celli)
    {
        const scalar csep = CSEP_.primitiveField()[celli];
        Cpsi.primitiveFieldRef()[celli] = scalar(0.2454)*std::pow(csep, scalar(-0.803));
    }

    forAll(Cpsi.boundaryField(), patchi)
    {
        const fvPatchScalarField& CSEPp = CSEP_.boundaryField()[patchi];
        fvPatchScalarField& Cpsip = Cpsi.boundaryFieldRef()[patchi];

        forAll(Cpsip, facei)
        {
            Cpsip[facei] = scalar(0.2454)*std::pow(CSEPp[facei], scalar(-0.803));
        }
    }

    // psi = C_psi (CSEP-1) Fblend            (A9)
    const volScalarField psi(Cpsi*(CSEP_ - scalar(1))*FBlend);

    // F_SEP = 1 + f_D * psi                  (A9)
    const volScalarField FSEP(scalar(1) + fD*psi);

    // Jet function F_JET                      (A14)
    volScalarField minSO(min(S, OmegaMag));
    minSO = min(minSO, omega_);

    const volScalarField Soo(minSO/(scalar(0.3)*omega_));

    const volScalarField xJet
    (
        scalar(4.0)
       *(
            (Soo - scalar(0.9))/scalar(0.15)
          - scalar(0.5)
        )
    );

    const volScalarField FJET
    (
        scalar(0.5)*(scalar(1) + tanh(xJet))
    );

    // F_MIX                                   (A13)
    const volScalarField FMIX
    (
        FSEP
      + (
            CMIX_
          + scalar(0.13)*CJET_*(FJET - scalar(1))
        )*(scalar(1) - FBlend)
    );


    // Cross diffusion CD                       (A6)
    const volScalarField CD
    (
        rho*(scalar(2)/sigmaOmega)
       *(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );

    // Destruction proxy                         (A12)
    const volScalarField Destr(Comega2_*FMIX*rho*sqr(omega_));

    // F_Lim                                    (A11)
    const volScalarField FLim(Destr/max(Destr, scalar(2)*CD));

    // F_NW                                     (A11), piecewise
    volScalarField FNW(scalar(0)*FBlend);

    // internal
    forAll(FNW.primitiveField(), celli)
    {
        const scalar cd = CD.primitiveField()[celli];

        if (cd > scalar(0))
        {
            FNW.primitiveFieldRef()[celli] = scalar(1);
        }
        else
        {
            const scalar cnw = CNW_.primitiveField()[celli];
            const scalar fd = fD.primitiveField()[celli];
            const scalar flim = FLim.primitiveField()[celli];

            FNW.primitiveFieldRef()[celli] =
                (
                    cnw
                  - (Cc1_.value() + Cc2_.value()*cnw)*fd
                )*flim;
        }
    }

    // patches
    forAll(FNW.boundaryField(), patchi)
    {
        const fvPatchScalarField& CDp = CD.boundaryField()[patchi];
        const fvPatchScalarField& CNWp = CNW_.boundaryField()[patchi];
        const fvPatchScalarField& fDp = fD.boundaryField()[patchi];
        const fvPatchScalarField& FLimp = FLim.boundaryField()[patchi];

        fvPatchScalarField& FNWp = FNW.boundaryFieldRef()[patchi];

        forAll(FNWp, facei)
        {
            const scalar cd = CDp[facei];

            if (cd > scalar(0))
            {
                FNWp[facei] = scalar(1);
            }
            else
            {
                const scalar cnw = CNWp[facei];
                const scalar fd = fDp[facei];
                const scalar flim = FLimp[facei];

                FNWp[facei] =
                    (
                        cnw
                      - (Cc1_.value() + Cc2_.value()*cnw)*fd
                    )*flim;
            }
        }
    }


    // ---------------------------------------------------------------------
    // Production terms
    // ---------------------------------------------------------------------

    // Raw production per-mass: G = nu_t * S^2 (A5), or nu_t*S*Omega (Eq.15)
    volScalarField G
    (
        katoLaunder_
      ? (this->nut_*S*OmegaMag)
      : (this->nut_*S2)
    );

    // Limited production per-mass: G~ = min(G, CPKlim*Cmu*k*omega) (A15/Eq.16)
    volScalarField Gtilde(G);

    if (productionLimiter_)
    {
        Gtilde = min(Gtilde, CPKlim_*Cmu_*k_*omega_);
    }

    // Omega production per-mass (A4):
    //   Pomega/rho = G / (nu_t + 0.001*nu)
    const volScalarField Pomega_m
    (
        G/(this->nut_ + scalar(0.001)*this->nu())
    );


    // ---------------------------------------------------------------------
    // k-equation (A1)
    // ---------------------------------------------------------------------

    const volScalarField DkEff(this->nu() + this->nut_/sigmaK);

    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff, k_)
     ==
        alpha()*rho()*Gtilde
      - fvm::Sp(alpha()*rho()*Cmu_*omega_(), k_)           // -Cmu*rho*k*omega
      + alpha()*rho()*Ck_*k_*omega_                        // +Ck*rho*k*omega (ML note)
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    kEqn.ref().solve();
    fvOptions.correct(k_);
    bound(k_, this->kMin_);


    // ---------------------------------------------------------------------
    // omega-equation (A2)
    // ---------------------------------------------------------------------

    const volScalarField DomegaEff(this->nu() + this->nut_/sigmaOmega);

    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(alpha, rho, omega_)
      + fvm::div(alphaRhoPhi, omega_)
      - fvm::laplacian(alpha*rho*DomegaEff, omega_)
     ==
        alpha()*rho()*Comega1*Pomega_m
      - fvm::Sp(alpha()*rho()*Comega2_*FMIX*omega_(), omega_) // -Cw2*FMIX*rho*omega^2
      + alpha()*FNW*CD                                        // +FNW*CD
      + alpha()*rho()*Comega_*sqr(omega_)                      // +Comega*rho*omega^2 (ML note)
      + fvOptions(alpha, rho, omega_)
    );

    omegaEqn.ref().relax();

    // Keep exact SST omega wall-function compatibility
    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());

    fvOptions.constrain(omegaEqn.ref());
    omegaEqn.ref().solve();
    fvOptions.correct(omega_);
    bound(omega_, this->omegaMin_);


    // Update nut with realizability limiter (A16)
    correctNut(S);
}

} // End namespace Foam

// ************************************************************************* //
