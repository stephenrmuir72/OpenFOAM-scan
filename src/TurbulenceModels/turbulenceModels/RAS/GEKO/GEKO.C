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
    const volSymmTensorField S(symm(fvc::grad(this->U_)));
    return tmp<volScalarField>::New("S2", 2.0*magSqr(S));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> GEKO<BasicTurbulenceModel>::Smag() const
{
    return tmp<volScalarField>::New("Smag", Foam::sqrt(max(S2(), VSMALL)));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> GEKO<BasicTurbulenceModel>::Omag() const
{
    const volTensorField W(skew(fvc::grad(this->U_)));
    return tmp<volScalarField>::New("Omag", sqrt(2.0*magSqr(W) + VSMALL));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> GEKO<BasicTurbulenceModel>::fD() const
{
    // (MM:Eq. A10)
    volScalarField ytilde
    (
        (1.0/kappa_)*(k_/max(omega_, this->omegaMin_))/this->nu()
    );

    return tmp<volScalarField>::New
    (
        "fD",
        1.0/(1.0 + sqr(ytilde/A1_))
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> GEKO<BasicTurbulenceModel>::FBlend() const
{
    if (wallDistanceFree_)
    {
        return tmp<volScalarField>::New("FBlend", scalar(1.0)*k_);
    }

    volScalarField kTilde("kTilde", max(k_, CFbLam_*this->nu()*omega_));
    volScalarField LT("LT", kTilde/(Cmu_*max(omega_, this->omegaMin_)));

    volScalarField xblend
    (
        "xblend",
        CFbTurb_*LT/max(y_, dimensionedScalar(dimLength,SMALL))
    );

    // (MM:Eq. A8)
    return tmp<volScalarField>::New("FBlend", tanh(pow4(xblend)));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> GEKO<BasicTurbulenceModel>::FSEP
(
    const volScalarField& Fblend
) const
{
    tmp<volScalarField> tCpsi = 0.2454*pow(CSEP_, -0.803);

    tmp<volScalarField> tpsi = tCpsi*(CSEP_ - 1.0)*Fblend;

    // (MM:Eq. A9)
    return tmp<volScalarField>::New("FSEP", 1.0 - fD()*tpsi);
}


template<class BasicTurbulenceModel>
tmp<volScalarField> GEKO<BasicTurbulenceModel>::FJET
(
    const volScalarField& Sm,
    const volScalarField& Om
) const
{
    tmp<volScalarField> tSoomin =
        min(min(Sm, Om), omega_)/(0.3*max(omega_, this->omegaMin_));

    tmp<volScalarField> txjet =
        (4.0*tSoomin() - 0.9)/(0.15 - 0.5*tSoomin() + SMALL);

    // (MM:Eq. A14)
    return tmp<volScalarField>::New("FJET", 0.5*(1.0 + tanh(txjet)));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> GEKO<BasicTurbulenceModel>::FMIX
(
    const volScalarField& Fsep,
    const volScalarField& Fblend,
    const volScalarField& Fjet
) const
{
    // (MM:Eq. A13)
    return tmp<volScalarField>::New
    (
        "FMIX",
        Fsep + (CMIX_ + 0.13*CJET_*Fjet - 1.0)*(1.0 - Fblend)
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
    tmp<volScalarField> tFLim = Destr/max(Destr, 2.0*CD);

    tmp<volScalarField> tfnw = (CNW_ - Cc1_ + Cc2_*CNW_)*Fd*tFLim;
    volScalarField& fnw = tfnw.ref();

    // (MM:Eq. A11)
    forAll (fnw, i)
    {
        if (CD[i] > 0)
        {
            fnw[i] = 1.0;
        }
    }
    fnw.correctBoundaryConditions();

    return tfnw;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> GEKO<BasicTurbulenceModel>::sigmaK
(
    const volScalarField& Fblend
) const
{
    // (MM:Eq. A17)
    return tmp<volScalarField>::New
    (
        "sigmaK",
        sigmaKtilde_*(1.0 + (CSEP_ - 1.0)*0.25*Fblend)
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> GEKO<BasicTurbulenceModel>::sigmaOmega
(
    const volScalarField& Fblend
) const
{
    // (MM:Eq. A17)
    return tmp<volScalarField>::New
    (
        "sigmaOmega",
        sigmaOmegaTilde_*(1.0 + (CSEP_ - 1.0)*Fblend)
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> GEKO<BasicTurbulenceModel>::Cw1
(
    const volScalarField& sigmaOmega
) const
{
    // (MM:Eq. A18)
    return tmp<volScalarField>::New
    (
        "Cw1",
        (1.0/Cmu_)*(Cw2_ - sqr(kappa_)*sqrt(Cmu_)/sigmaOmega)
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> GEKO<BasicTurbulenceModel>::CDterm
(
    const volScalarField& sigmaOmega
) const
{
    return tmp<volScalarField>::New
    (
        "CD",
        0.5/sigmaOmega
       *((fvc::grad(k_) & fvc::grad(omega_))/max(omega_, this->omegaMin_))
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> GEKO<BasicTurbulenceModel>::limitedPk
(
    const volScalarField& Pk
) const
{
    tmp<volScalarField> tepsilon = Cmu_*k_*omega_;

    // (MM: A15)
    return tmp<volScalarField>::New("PkLimited", min(Pk, CPKlim_*tepsilon));
}


template<class BasicTurbulenceModel>
void GEKO<BasicTurbulenceModel>::realizabilityLimitNut(volScalarField& nut)
{
    nut = min(nut, CRealize_*k_/max(Smag(), VSMALL));
}


template<class BasicTurbulenceModel>
void GEKO<BasicTurbulenceModel>::correctNut()
{
    // (MM:Eq. A16)
    this->nut_ = k_/max(omega_, this->omegaMin_);

    realizabilityLimitNut(this->nut_);

    // bound(this->nut_, 0.0);  PENDING
    this->nut_.correctBoundaryConditions();

    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


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

    // Construct local convenience references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    const volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();


    //
    const volScalarField Sm(Smag());
    const volScalarField Om(Omag());


    // Calculate GEKO functions
    tmp<volScalarField> tFblend = FBlend();
    const volScalarField& Fblend = tFblend();

    tmp<volScalarField> tFsep = FSEP(Fblend);
    const volScalarField& Fsep = tFsep();

    tmp<volScalarField> tFjet = FJET(Sm, Om);
    const volScalarField& Fjet = tFjet();

    tmp<volScalarField> tFmix = FMIX(Fsep, Fblend, Fjet);
    const volScalarField& Fmix = tFmix();

    tmp<volScalarField> tsigK = sigmaK(Fblend);
    const volScalarField& sigK = tsigK();

    tmp<volScalarField> tsigO = sigmaOmega(Fblend);
    const volScalarField& sigO = tsigO();

    tmp<volScalarField> tCw1 = Cw1(sigO);
    const volScalarField& Cw1f = tCw1();

    tmp<volScalarField> tCD = CDterm(sigO);
    const volScalarField& CD = tCD();

    volScalarField Destr("Destr", Cw2_*Fmix*sqr(omega_));

    tmp<volScalarField> tFNW = FNW(CD, fD(), Destr);
    const volScalarField& FNWf = tFNW();


    // Apply the realizability limiter onto nut
    realizabilityLimitNut(this->nut_);

    tmp<volScalarField> tPk;
    if (useKatoLaunder_)
    {
        tPk = tmp<volScalarField>::New("Pk", this->nut_*Sm*Om);
    }
    else
    {
        tPk = tmp<volScalarField>::New("Pk", this->nut_*S2());
    }
    const volScalarField& Pk = tPk.cref();


    // Production limiter (A15)
    tmp<volScalarField> tPkLim = limitedPk(Pk);
    const volScalarField& PkLim = tPkLim();

    volScalarField Pomega
    (
        "Pomega",
        (this->nut_ + 0.001*this->nu())*PkLim/max(k_, this->kMin_)
    );


    tmp<volTensorField> tgradU = fvc::grad(U);
    const volScalarField::Internal GbyNu
    (
        tgradU().v() && devTwoSymm(tgradU().v())
    );
    const volScalarField::Internal G(this->GName(), nut()*GbyNu);
    tgradU.clear();

    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();
    // Push any changed cell values to coupled neighbours
    omega_.boundaryFieldRef().template evaluateCoupled<coupledFvPatch>();


    // Turbulence specific dissipation rate equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(alpha, rho, omega_)
      + fvm::div(alphaRhoPhi, omega_)
      - fvm::laplacian(alpha*rho*DomegaEff(), omega_)
     ==
        alpha()*rho()*Cw1f*Pomega
      + alpha()*rho()*FNWf*CD
      - fvm::Sp(Cw2_*alpha()*rho()*Fmix()*omega_(), omega_)
      + fvOptions(alpha, rho, omega_)
    );

    omegaEqn.ref().relax();
    fvOptions.constrain(omegaEqn.ref());
    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());  // PENDING?
    solve(omegaEqn);
    fvOptions.correct(omega_);
    bound(omega_, this->omegaMin_);


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

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);


    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
