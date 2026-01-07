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

\*---------------------------------------------------------------------------*/

#include "GEKO.H"
#include "fvOptions.H"
#include "bound.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> GEKO<BasicTurbulenceModel>::S2() const
{
    DebugVar("S2-1");
    const volVectorField& U = this->U_;
    tmp<volTensorField> tgradU = fvc::grad(U);  // [1/s]

    return tmp<volScalarField>::New
    (
        2.0*magSqr(symm(tgradU))  // [1/s^2]
    );
    DebugVar("S2-2");
}


template<class BasicTurbulenceModel>
tmp<volScalarField> GEKO<BasicTurbulenceModel>::Smag() const
{
    DebugVar("Smag-1");
    const dimensionedScalar eps(dimless/sqr(dimTime), VSMALL);
    return tmp<volScalarField>::New(Foam::sqrt(max(S2(), eps)));  // [1/s]
}


template<class BasicTurbulenceModel>
tmp<volScalarField> GEKO<BasicTurbulenceModel>::Omag() const
{
    DebugVar("Omag-1");
    const volVectorField& U = this->U_;
    tmp<volTensorField> tgradU = fvc::grad(U);  // [1/s]
    const dimensionedScalar eps(dimless/sqr(dimTime), VSMALL);

    return tmp<volScalarField>::New
    (
        sqrt(2.0*magSqr(skew(tgradU)) + eps)  // [1/s]
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> GEKO<BasicTurbulenceModel>::fD() const
{
    DebugVar("fD-1");
    // (MM:Eq. A10)
    volScalarField ytilde  // [dimless]
    (
        (1.0/kappa_)*(k_/max(omega_, this->omegaMin_))/this->nu()
    );
    DebugVar("fD-2");

    return tmp<volScalarField>::New
    (
        1.0/(1.0 + sqr(ytilde/A1_))  // [dimless]
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> GEKO<BasicTurbulenceModel>::FBlend() const
{
    DebugVar("FBlend-1");
    if (wallDistanceFree_)
    {
    DebugVar("FBlend-2");
        return tmp<volScalarField>::New(scalar(1.0)*k_);  // [m^2/s^2]
    }

    tmp<volScalarField> tkTilde =  // [m^2/s^2]
        max(k_, CFbLam_*this->nu()*omega_);

    tmp<volScalarField> tLT =  // [m]
        Foam::sqrt(tkTilde)/(Cmu_*max(omega_, this->omegaMin_));

    tmp<volScalarField> txBlend =  // [dimless]
        CFbTurb_*tLT/max(y_, dimensionedScalar(dimLength, SMALL));

    DebugVar("FBlend-3");

    // (MM:Eq. A8)
    return tmp<volScalarField>::New(tanh(pow4(txBlend)));  // [dimless]
}


template<class BasicTurbulenceModel>
tmp<volScalarField> GEKO<BasicTurbulenceModel>::FSEP
(
    const volScalarField& Fblend
) const
{
    DebugVar("FSEP-1");
    tmp<volScalarField> tCpsi = 0.2454*pow(CSEP_, -0.803);

    tmp<volScalarField> tpsi = tCpsi*(CSEP_ - 1.0)*Fblend;  // [dimless]
    DebugVar("FSEP-2");

    // (MM:Eq. A9)
    return tmp<volScalarField>::New(1.0 - fD()*tpsi);  // [dimless]
}


template<class BasicTurbulenceModel>
tmp<volScalarField> GEKO<BasicTurbulenceModel>::FJET
(
    const volScalarField& Sm,
    const volScalarField& Om
) const
{
    DebugVar("FJET-1");
    tmp<volScalarField> tSoomin =  // [dimless]
        min(min(Sm, Om), omega_)/(0.3*max(omega_, this->omegaMin_));

    const dimensionedScalar eps(dimless, SMALL);

    DebugVar("FJET-2");
    tmp<volScalarField> txjet =  // [dimless]
        (4.0*tSoomin() - 0.9)/(0.15 - 0.5*tSoomin() + eps);

    DebugVar("FJET-3");
    // (MM:Eq. A14)
    return tmp<volScalarField>::New(0.5*(1.0 + tanh(txjet)));  // [dimless]
}


template<class BasicTurbulenceModel>
tmp<volScalarField> GEKO<BasicTurbulenceModel>::FMIX
(
    const volScalarField& Fsep,
    const volScalarField& Fblend,
    const volScalarField& Fjet
) const
{
    DebugVar("FMIX-1");
    // (MM:Eq. A13)
    return tmp<volScalarField>::New
    (
        Fsep + (CMIX_ + 0.13*CJET_*Fjet - 1.0)*(1.0 - Fblend)  // [dimless]
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> GEKO<BasicTurbulenceModel>::FNW
(
    const volScalarField& CD,
    const volScalarField& Fd,
    const volScalarField& Destr
) const
{
    DebugVar("FNW-1");
    tmp<volScalarField> tFLim = Destr/max(Destr, 2.0*CD);  // [1/s^2]

    tmp<volScalarField> tfnw = (CNW_ - Cc1_ + Cc2_*CNW_)*Fd*tFLim;  // [1/s^2]
    volScalarField& fnw = tfnw.ref();
    DebugVar("FNW-2");

    // const dimensionedScalar one(fnw.dimensions(), 1);

    // (MM:Eq. A11)
    forAll (fnw, i)
    {
        if (CD[i] > 0)
        {
//            fnw[i] = one;
            fnw[i] = 1.0;
        }
    }
    fnw.correctBoundaryConditions();
    DebugVar("FNW-3");

    return tfnw;  // [1/s^2]
}


template<class BasicTurbulenceModel>
tmp<volScalarField> GEKO<BasicTurbulenceModel>::sigmaK
(
    const volScalarField& Fblend
) const
{
    DebugVar("sigmaK-1");
    // (MM:Eq. A17)
    return tmp<volScalarField>::New  // [dimless]
    (
        sigmaKtilde_*(1.0 + (CSEP_ - 1.0)*0.25*Fblend)
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> GEKO<BasicTurbulenceModel>::sigmaOmega
(
    const volScalarField& Fblend
) const
{
    DebugVar("sigmaOmega-1");
    // (MM:Eq. A17)
    return tmp<volScalarField>::New  // [dimless]
    (
        sigmaOmegaTilde_*(1.0 + (CSEP_ - 1.0)*Fblend)
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> GEKO<BasicTurbulenceModel>::Cw1
(
    const volScalarField& sigmaOmega
) const
{
    DebugVar("Cw1 - 1");
    // (MM:Eq. A18)
    return tmp<volScalarField>::New  // [dimless]
    (
        (1.0/Cmu_)*(Cw2_ - sqr(kappa_)*sqrt(Cmu_)/sigmaOmega)
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> GEKO<BasicTurbulenceModel>::CDterm
(
    const volScalarField& sigmaOmega
) const
{
    DebugVar("CDterm-1");
    // (MM:Eq. A6)
    return tmp<volScalarField>::New  // [1/s^2]
    (
        0.5/sigmaOmega
       *(
        (fvc::grad(k_) & fvc::grad(omega_))
       /max(omega_, this->omegaMin_)
        )
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> GEKO<BasicTurbulenceModel>::limitedPk
(
    const volScalarField& Pk
) const
{
    DebugVar("limitedPk-1");
    tmp<volScalarField> tepsilon = Cmu_*k_*omega_;  // [m^2/s^3]

    // (MM: A15)
    return tmp<volScalarField>::New(min(Pk, CPKlim_*tepsilon));  // [m^2/s^3]
}


template<class BasicTurbulenceModel>
void GEKO<BasicTurbulenceModel>::realizabilityLimitNut(volScalarField& nut)
{
    DebugVar("realizability-1");
    tmp<volScalarField> tSmag = Smag();

    const dimensionedScalar eps(tSmag().dimensions(), VSMALL);

    // (MM:Eq. A16)
    nut = min(nut, CRealize_*k_/max(tSmag, eps));  // [m^2/s]
    DebugVar("realizability-2");
}


template<class BasicTurbulenceModel>
void GEKO<BasicTurbulenceModel>::correctNut()
{
    DebugVar("correctNut-1");
    // (MM:Eq. A16)
    this->nut_ = k_/max(omega_, this->omegaMin_);  // [m^2/s]

    realizabilityLimitNut(this->nut_);

    // bound(this->nut_, 0.0);  PENDING
    this->nut_.correctBoundaryConditions();

    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
    DebugVar("correctNut-2");
}


/*
template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> GEKO<BasicTurbulenceModel>::kSource() const
{
    return tmp<fvScalarMatrix>::New
    (
        k_,
        dimVolume*this->rho_.dimensions()*k_.dimensions()/dimTime
    );
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> GEKO<BasicTurbulenceModel>::omegaSource() const
{
    return tmp<fvScalarMatrix>::New
    (
        omega_,
        dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
    );
}
*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
GEKO<BasicTurbulenceModel>::GEKO
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel>>
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

    y_(wallDist::New(this->mesh_).y()),
    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
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
            IOobject::groupName("omega", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    CSEP_
    (
        IOobject
        (
            IOobject::groupName("CSEP", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimless, 1.75)
    ),
    CNW_
    (
        IOobject
        (
            IOobject::groupName("CNW", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimless, 0.0)
    ),
    CMIX_
    (
        IOobject
        (
            IOobject::groupName("CMIX", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimless, 0.5)
    ),
    CJET_
    (
        IOobject
        (
            IOobject::groupName("CJET", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimless, 1.0)
    ),

    Cmu_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    Cw2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cw2",
            this->coeffDict_,
            0.083
        )
    ),
    sigmaKtilde_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmaKtilde",
            this->coeffDict_,
            1.0
        )
    ),
    sigmaOmegaTilde_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmaOmegaTilde",
            this->coeffDict_,
            1.17
        )
    ),
    kappa_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "kappa",
            this->coeffDict_,
            0.41
        )
    ),
    CFbTurb_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "CFbTurb",
            this->coeffDict_,
            2.0
        )
    ),
    CFbLam_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "CFbLam",
            this->coeffDict_,
            1.0
        )
    ),
    A1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A1",
            this->coeffDict_,
            15.0
        )
    ),
    Cc1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cc1",
            this->coeffDict_,
            1.7
        )
    ),
    Cc2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cc2",
            this->coeffDict_,
            1.4
        )
    ),

    wallDistanceFree_
    (
        Switch::getOrAddToDict
        (
            "wallDistanceFree",
            this->coeffDict_,
            false
        )
    ),
    useKatoLaunder_
    (
        Switch::getOrAddToDict
        (
            "useKatoLaunder",
            this->coeffDict_,
            true
        )
    ),
    CPKlim_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "CPKlim",
            this->coeffDict_,
            10.0
        )
    ),
    CRealize_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "CRealize",
            this->coeffDict_,
            0.577
        )
    )
{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool GEKO<BasicTurbulenceModel>::read()
{
    if (!eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        return false;
    }

    Cmu_.readIfPresent(this->coeffDict());
    Cw2_.readIfPresent(this->coeffDict());
    sigmaKtilde_.readIfPresent(this->coeffDict());
    sigmaOmegaTilde_.readIfPresent(this->coeffDict());
    kappa_.readIfPresent(this->coeffDict());
    CFbTurb_.readIfPresent(this->coeffDict());
    CFbLam_.readIfPresent(this->coeffDict());
    A1_.readIfPresent(this->coeffDict());
    Cc1_.readIfPresent(this->coeffDict());
    Cc2_.readIfPresent(this->coeffDict());

    wallDistanceFree_.readIfPresent("wallDistanceFree", this->coeffDict());
    useKatoLaunder_.readIfPresent("useKatoLaunder", this->coeffDict());
    CPKlim_.readIfPresent(this->coeffDict());
    CRealize_.readIfPresent(this->coeffDict());

    return true;
}


template<class BasicTurbulenceModel>
void GEKO<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_) return;

    DebugVar("loop-1");
    // Construct local convenience references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    const volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    DebugVar("loop-2");

    //
    const volScalarField Sm(Smag());
    const volScalarField Om(Omag());

    DebugVar("loop-3");

    // Calculate GEKO functions
    tmp<volScalarField> tFblend = FBlend();
    const volScalarField& Fblend = tFblend();
    DebugVar("loop-4");

    tmp<volScalarField> tFsep = FSEP(Fblend);
    const volScalarField& Fsep = tFsep();
    DebugVar("loop-5");

    tmp<volScalarField> tFjet = FJET(Sm, Om);
    const volScalarField& Fjet = tFjet();
    DebugVar("loop-6");

    tmp<volScalarField> tFmix = FMIX(Fsep, Fblend, Fjet);
    const volScalarField& Fmix = tFmix();
    DebugVar("loop-7");

    tmp<volScalarField> tsigK = sigmaK(Fblend);
    const volScalarField& sigK = tsigK();
    DebugVar("loop-8");

    tmp<volScalarField> tsigO = sigmaOmega(Fblend);
    const volScalarField& sigO = tsigO();
    DebugVar("loop-9");

    tmp<volScalarField> tCw1 = Cw1(sigO);
    const volScalarField& Cw1f = tCw1();
    DebugVar("loop-10");

    tmp<volScalarField> tCD = CDterm(sigO);
    const volScalarField& CD = tCD();
    DebugVar("loop-11");

    const volScalarField Destr(Cw2_*Fmix*sqr(omega_));  // [1/s^2]
    DebugVar("loop-12");

    tmp<volScalarField> tFNW = FNW(CD, fD(), Destr);
    const volScalarField& FNWf = tFNW();
    DebugVar("loop-13");


    // Apply the realizability limiter onto nut
    //realizabilityLimitNut(this->nut_);
    DebugVar("loop-14");

    tmp<volScalarField> tPk;
    if (useKatoLaunder_)
    {
        tPk = tmp<volScalarField>::New(this->nut_*Sm*Om);
    }
    else
    {
        tPk = tmp<volScalarField>::New(this->nut_*S2());
    }
    const volScalarField& Pk = tPk.cref();  // [m^2/s^3]
    DebugVar("loop-15");


    // Production limiter (A15)
    tmp<volScalarField> tPkLim = limitedPk(Pk);
    const volScalarField& PkLim = tPkLim();
    DebugVar("loop-16");

    volScalarField Pomega  // [1/s^2]
    (
        // (this->nut_ + 0.001*this->nu())*PkLim/max(k_, this->kMin_)
        PkLim/(this->nut_ + 0.001*this->nu())
    );

    DebugVar("loop-16");

    tmp<volTensorField> tgradU = fvc::grad(U);
    const volScalarField::Internal GbyNu
    (
        tgradU().v() && devTwoSymm(tgradU().v())
    );
    const volScalarField::Internal G(this->GName(), nut()*GbyNu);
    DebugVar("loop-17");

    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();
    // Push any changed cell values to coupled neighbours
    omega_.boundaryFieldRef().template evaluateCoupled<coupledFvPatch>();

    DebugVar("loop-18");


    // Turbulence specific dissipation rate equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(alpha, rho, omega_)   // [1/s^2]
      + fvm::div(alphaRhoPhi, omega_)  // [1/s^2]
      - fvm::laplacian(alpha*rho*DomegaEff(), omega_)  // [1/s^2]
     ==
        alpha()*rho()*Cw1f*Pomega  // [1/s^2]
      + alpha()*rho()*FNWf*CD
      - fvm::Sp(Cw2_*alpha()*rho()*Fmix()*omega_(), omega_)
      + fvOptions(alpha, rho, omega_)
    );
    DebugVar("loop-19");

    omegaEqn.ref().relax();
    DebugVar("loop-20");
    fvOptions.constrain(omegaEqn.ref());
    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());  // PENDING?
    DebugVar("loop-21");
    solve(omegaEqn);
    fvOptions.correct(omega_);
    DebugVar("loop-22");
    bound(omega_, this->omegaMin_);
    DebugVar("loop-23");


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha()*rho()*PkLim
      // - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
      - fvm::Sp(Cmu_*alpha()*rho()*omega_(), k_)
      + fvOptions(alpha, rho, k_)
    );
    DebugVar("loop-24");

    kEqn.ref().relax();
    DebugVar("loop-25");
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    DebugVar("loop-26");
    fvOptions.correct(k_);
    bound(k_, this->kMin_);
    DebugVar("loop-27");


    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
