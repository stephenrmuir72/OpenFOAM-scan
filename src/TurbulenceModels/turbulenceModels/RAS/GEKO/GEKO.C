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
void GEKO<BasicTurbulenceModel>::correctNut()
{

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


    wallDistanceFree_
    (
        Switch::getOrAddToDict
        (
            "wallDistanceFree",
            this->coeffDict_,
            false
        )
    ),
    productionLimiter_
    (
        Switch::getOrAddToDict
        (
            "productionLimiter",
            this->coeffDict_,
            true
        )
    ),
    katoLaunder_
    (
        Switch::getOrAddToDict
        (
            "katoLaunder",
            this->coeffDict_,
            false
        )
    ),
    dilatationCorrection_
    (
        Switch::getOrAddToDict
        (
            "dilatationCorrection",
            this->coeffDict_,
            false
        )
    ),


    Cmu_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    Comega2_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Comega2",
            this->coeffDict_,
            0.083
        )
    ),
    sigmakTilde_
    (
        dimensionedScalar::getOrAddToDict
        (
            "sigmakTilde",
            this->coeffDict_,
            1.0
        )
    ),
    sigmaOmegaTilde_
    (
        dimensionedScalar::getOrAddToDict
        (
            "sigmaOmegaTilde",
            this->coeffDict_,
            1.17
        )
    ),
    CRealize_
    (
        dimensionedScalar::getOrAddToDict
        (
            "CRealize",
            this->coeffDict_,
            scalar(1.0)/std::sqrt(scalar(3.0))  // (MM:Eq. 17)
        )
    ),
    kappa_
    (
        dimensionedScalar::getOrAddToDict
        (
            "kappa",
            this->coeffDict_,
            0.41
        )
    ),
    CFbTurb_
    (
        dimensionedScalar::getOrAddToDict
        (
            "CFbTurb",
            this->coeffDict_,
            2.0
        )
    ),
    CFbLam_
    (
        dimensionedScalar::getOrAddToDict
        (
            "CFbLam",
            this->coeffDict_,
            1.0  // increase to 25 in case of transition simulations (MM:p. 4)
        )
    ),
    Aplus_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Aplus",
            this->coeffDict_,
            15.0
        )
    ),
    Cc1_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Cc1",
            this->coeffDict_,
            1.7
        )
    ),
    Cc2_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Cc2",
            this->coeffDict_,
            1.4
        )
    ),
    CPKlim_
    (
        dimensionedScalar::getOrAddToDict
        (
            "CPKlim",
            this->coeffDict_,
            10.0
        )
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
        dimensionedScalar(dimless, 0.5)
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
        dimensionedScalar(dimless, 0.0)
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


    Ck_
    (
        IOobject
        (
            IOobject::groupName("Ck", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimless, 0.0)
    ),
    Comega_
    (
        IOobject
        (
            IOobject::groupName("Comega", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimless, 0.0)
    ),


    FSEP_
    (
        IOobject
        (
            IOobject::groupName("FSEP", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimless, 1.0)
    ),
    FNW_
    (
        IOobject
        (
            IOobject::groupName("FNW", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimless, 1.0)
    ),
    FMIX_
    (
        IOobject
        (
            IOobject::groupName("FMIX", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimless, 1.0)
    ),
    FJET_
    (
        IOobject
        (
            IOobject::groupName("FJET", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimless, 1.0)
    ),
    Fblend_
    (
        IOobject
        (
            IOobject::groupName("Fblend", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimless, 1.0)  // wall-distance free variant
    )
{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }

    if (Cmu_.value() < SMALL)
    {
        FatalIOErrorInFunction(this->coeffDict_)
            << "Invalid value for Cmu: " << Cmu_.value()
            << ". Valid range is (0,inf]."
            << exit(FatalIOError);
    }

    if (kappa_.value() < SMALL)
    {
        FatalIOErrorInFunction(this->coeffDict_)
            << "Invalid value for kappa: " << kappa_.value()
            << ". Valid range is (0,inf]."
            << exit(FatalIOError);
    }

    if (Aplus_.value() < SMALL)
    {
        FatalIOErrorInFunction(this->coeffDict_)
            << "Invalid value for Aplus: " << Aplus_.value()
            << ". Valid range is (0,inf]."
            << exit(FatalIOError);
    }

    // y_ cannot go to zero or negative.

    // omega_ and nu_ should not go to zero anytime.
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool GEKO<BasicTurbulenceModel>::read()
{
    if (!eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        return false;
    }

    wallDistanceFree_.readIfPresent("wallDistanceFree", this->coeffDict());
    productionLimiter_.readIfPresent("productionLimiter", this->coeffDict());
    katoLaunder_.readIfPresent("katoLaunder", this->coeffDict());
    dilatationCorrection_.readIfPresent
    (
        "dilatationCorrection",
        this->coeffDict()
    );

    Cmu_.readIfPresent(this->coeffDict());
    Comega2_.readIfPresent(this->coeffDict());
    sigmakTilde_.readIfPresent(this->coeffDict());
    sigmaOmegaTilde_.readIfPresent(this->coeffDict());
    CRealize_.readIfPresent(this->coeffDict());
    kappa_.readIfPresent(this->coeffDict());
    CFbTurb_.readIfPresent(this->coeffDict());
    CFbLam_.readIfPresent(this->coeffDict());
    Aplus_.readIfPresent(this->coeffDict());
    Cc1_.readIfPresent(this->coeffDict());
    Cc2_.readIfPresent(this->coeffDict());
    CPKlim_.readIfPresent(this->coeffDict());

    if (Cmu_.value() < SMALL)
    {
        FatalIOErrorInFunction(this->coeffDict_)
            << "Invalid value for Cmu: " << Cmu_.value()
            << ". Valid range is (0,inf]."
            << exit(FatalIOError);
    }

    if (kappa_.value() < SMALL)
    {
        FatalIOErrorInFunction(this->coeffDict_)
            << "Invalid value for kappa: " << kappa_.value()
            << ". Valid range is (0,inf]."
            << exit(FatalIOError);
    }

    if (Aplus_.value() < SMALL)
    {
        FatalIOErrorInFunction(this->coeffDict_)
            << "Invalid value for Aplus: " << Aplus_.value()
            << ". Valid range is (0,inf]."
            << exit(FatalIOError);
    }

    return true;
}


template<class BasicTurbulenceModel>
void GEKO<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Construct local convenience references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    const volScalarField& nut = this->nut_;

    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    // Calculate the velocity gradient tensor in Hessian form (delta_i u_j)
    // Tranpose of the classical Jacobian form (delta_j u_i)
    tmp<volTensorField> tgradU = fvc::grad(U);
    const volTensorField& gradU = tgradU.cref();

    // Strain-rate tensor (MM:Eq. 4)
    const volScalarField S2(2.0*magSqr(symm(gradU)));

    // Strain-rate tensor magnitude (MM:Eq. 4)
    const volScalarField S(sqrt(S2));

    // Rotation tensor magnitude (MM:Eq. A14)
    const volScalarField OmegaMag(sqrt(2.0*magSqr(skew(gradU)))));

    // Update omega and G at the wall  - PENDING G
    omega_.boundaryFieldRef().updateCoeffs();
    // Push any changed cell values to coupled neighbours
    omega_.boundaryFieldRef().template evaluateCoupled<coupledFvPatch>();


    // Calculate the blending function field (MM:Eq. A8)

    if (!wallDistanceFree_)
    {
        const volScalarField kTilde(max(k_, CFbLam_*this->nu()*omega_)); // PENDING: cannot go to zero or neg.

        const volScalarField LT(sqrt(kTilde)/(Cmu_*omega_));  // PENDING: omega cannot go to zero or neg.

        const volScalarField xBlend(CFbTurb_*LT/y_);  // PENDING: y_ cannot go to zero or neg.

        FBlend_ = tanh(pow4(xBlend));
        FBlend_.correctBoundaryConditions();
    }


    // Calculate the blending-function coefficient fields (MM:Eqs. A17-A18)

    const volScalarField sigmaK
    (
        sigmakTilde_*(scalar(1) + (CSEP_ - scalar(1))*scalar(0.25)*FBlend_)
    );

    const volScalarField sigmaOmega
    (
        sigmaOmegaTilde_*(scalar(1) + (CSEP_ - scalar(1))*FBlend_)  // PENDING: sigmaOmega cannot be zero
    );

    const scalar sqrtCmu = Foam::sqrt(Cmu_.value());
    const volScalarField Comega1
    (
        (scalar(1)/Cmu_)
       *(
            Comega2_
          - sqr(kappa_)*sqrtCmu/sigmaOmega
        )
    );


    // Calculate the flow-separation calibration field, FSEP (MM:Eqs. A9-10)

    const volScalarField yTildePlus
    (
        (scalar(1)/kappa_)*(k_/omega_)/this->nu()
    );

    const volScalarField fD
    (
        scalar(1)/(scalar(1) + sqr(yTildePlus/APlus_))
    );

    const volScalarField Cpsi
    (
        scalar(0.2454)*Foam::pow(CSEP_, scalar(-0.803))
    );

    const volScalarField psi(Cpsi*(CSEP_ - scalar(1))*FBlend_);

    FSEP_ = scalar(1) + fD*psi;
    FSEP_.correctBoundaryConditions();


    // Calculate the jet-flow calibration field, FJET (MM:Eqs. A14)

    const volScalarField Soo
    (
        min
        (
            min(S, OmegaMag),
            omega_
        )
       /(scalar(0.3)*omega_)
    );

    const volScalarField xJet
    (
        scalar(4.0)
       *(
            (Soo - scalar(0.9))/scalar(0.15)
          - scalar(0.5)
        )
    );

    FJET_ = scalar(0.5)*(scalar(1) + tanh(xJet));
    FJET_.correctBoundaryConditions();


    // Calculate the mixing function field, FMIX (MM:Eq. A13)

    FMIX_ =
        FSEP_
      + (
            CMIX_
          + scalar(0.13)*CJET_*(FJET_ - scalar(1))
        )
       *(scalar(1) - FBlend_);
    FMIX_.correctBoundaryConditions();


    // Calculate the non-equilibrium calibration field, FNW (MM:Eqs. A11)

    // Cross-diffusion term (MM:Eq. A6)  - PENDING rho??
    const volScalarField CD
    (
        rho*(scalar(2)/sigmaOmega)
       *(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );

    // Destruction term (MM:Eq. A12)
    const volScalarField Destr(Comega2_*FMIX*rho*sqr(omega_));

    const volScalarField FLim(Destr/max(Destr, scalar(2)*CD));

    forAll(FNW_, celli)
    {
        if (CD[celli] > scalar(0))
        {
            FNW_[celli] = scalar(1);
        }
        else
        {
            FNW_[celli] =
                CNW_[celli] - (Cc1_ + Cc2_*CNW_[celli])*fd[celli]*FLim[celli];
        }
    }
    FNW_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
