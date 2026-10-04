/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "liquidFilmModel.H"
#include "volFields.H"
#include "spatialSpreading.H"
#include "addToRunTimeSelectionTable.H"
#include "Map.H"

#include <cmath>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(liquidFilmModel, 0);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void liquidFilmModel::correctThermoFields()
{
    scalarField X(thermo_.size(), 1);

    forAll(rho_, faceI)
    {
        rho_[faceI] = thermo_.rho(pRef_, Tf_[faceI], X);
        mu_[faceI] = thermo_.mu(pRef_, Tf_[faceI], X);
        sigma_[faceI] = thermo_.sigma(pRef_, Tf_[faceI], X);
        Cp_[faceI] = thermo_.Cp(pRef_, Tf_[faceI], X);
    }

    forAll(regionMesh().boundary(), patchI)
    {
        const scalarField& patchTf = Tf_.boundaryFieldRef()[patchI];

        scalarField& patchRho = rho_.boundaryFieldRef()[patchI];
        scalarField& patchmu = mu_.boundaryFieldRef()[patchI];
        scalarField& patchsigma = sigma_.boundaryFieldRef()[patchI];
        scalarField& patchCp = Cp_.boundaryFieldRef()[patchI];

        forAll(patchRho, edgeI)
        {
            patchRho[edgeI] = thermo_.rho(pRef_, patchTf[edgeI], X);
            patchmu[edgeI] = thermo_.mu(pRef_, patchTf[edgeI], X);
            patchsigma[edgeI] = thermo_.sigma(pRef_, patchTf[edgeI], X);
            patchCp[edgeI] = thermo_.Cp(pRef_, patchTf[edgeI], X);
        }
    }

    //Initialize pf_
    pf_ = rho_*gn_*h_ - sigma_*fac::laplacian(h_);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

liquidFilmModel::liquidFilmModel
(
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    liquidFilmBase(modelType, mesh, dict),
    thermo_(dict.subDict("thermo")),
    rho_
    (
        IOobject
        (
            "rhof",
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimDensity, Zero)
    ),
    mu_
    (
        IOobject
        (
            "muf",
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimViscosity, Zero)
    ),
    Tf_
    (
        IOobject
        (
            "Tf_" + regionName_,
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimTemperature, Zero)
    ),
    Cp_
    (
        IOobject
        (
            "Cp_" + regionName_,
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimEnergy/dimTemperature, Zero)
    ),
    sigma_
    (
        IOobject
        (
            "sigmaf",
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimMass/sqr(dimTime), Zero)
    ),
    hRho_
    (
        IOobject
        (
            h_.name() + "*" + rho_.name(),
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(h_.dimensions()*rho_.dimensions(), Zero)
    ),
    rhoSp_
    (
        IOobject
        (
            "rhoSp",
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimVelocity, Zero)
    ),
    USp_
    (
        IOobject
        (
            "USp",
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedVector(sqr(dimVelocity), Zero)
    ),
    pnSp_
    (
        IOobject
        (
            "pnSp",
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimPressure, Zero)
    ),
    cloudMassTrans_
    (
        IOobject
        (
            "cloudMassTrans",
            primaryMesh().time().timeName(),
            primaryMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primaryMesh(),
        dimensionedScalar(dimMass, Zero)
    ),
    cloudDiameterTrans_
    (
        IOobject
        (
            "cloudDiameterTrans",
            primaryMesh().time().timeName(),
            primaryMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primaryMesh(),
        dimensionedScalar(dimLength, Zero)
    ),
    turbulence_(filmTurbulenceModel::New(*this, dict)),
    impacts_(),
    availableMass_(regionMesh().faces().size(), Zero),
    injection_(*this, dict),
    forces_(*this, dict)
{
    if (dict.readIfPresent("T0", Tref_))
    {
        Info<< "Initialised film temperature to T0" << endl;

        Tf_ = dimensionedScalar("T0", dimTemperature, dict);
    }
    correctThermoFields();

    spatialSpreadingPtr_ = spatialSpreading::NewIfPresent(regionMesh(), dict);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

liquidFilmModel::~liquidFilmModel() = default;


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const volScalarField& liquidFilmModel::cloudMassTrans() const
{
    return cloudMassTrans_;
}


const volScalarField& liquidFilmModel::cloudDiameterTrans() const
{
    return cloudDiameterTrans_;
}

void liquidFilmModel::addParcelSources
(
    const label patchi,
    const label facei,
    const scalar mass,
    const vector& momentum,
    const scalar normalImpulse,
    const scalar energy,
    const scalar diameter
)
{
    if (!std::isfinite(mass))
    {
        FatalErrorInFunction
            << "Invalid non-finite parcel mass on patch " << patchi
            << ", face " << facei << ": mass=" << mass
            << abort(FatalError);
    }

    if (mass > 0 && spatialSpreadingPtr_)
    {
        if
        (
            !std::isfinite(mag(momentum))
         || !std::isfinite(normalImpulse)
         || !(diameter > 0)
         || !std::isfinite(diameter)
        )
        {
            FatalErrorInFunction
                << "Invalid parcel impact on patch " << patchi
                << ", face " << facei << ": mass=" << mass
                << ", diameter=" << diameter
                << abort(FatalError);
        }

        impacts_.append
        (
            parcelImpact{patchi, facei, mass, momentum, normalImpulse, diameter}
        );
        // Spread this impact from its retained record, not the aggregate.
        return;
    }

    liquidFilmBase::addParcelSources
    (
        patchi, facei, mass, momentum, normalImpulse, energy, diameter
    );
}


void liquidFilmModel::preEvolveRegion()
{
    liquidFilmBase::preEvolveRegion();


    cloudMassTrans_ == Zero;
    cloudDiameterTrans_ == Zero;

    const scalar deltaT = primaryMesh().time().deltaTValue();
    const scalarField rAreaDeltaT(scalar(1)/deltaT/regionMesh().S().field());

    // Map the total mass, mom [kg.m/s] and pnSource from particles

    vsm().mapToSurface(massSource_, rhoSp_.primitiveFieldRef());

    vsm().mapToSurface(momentumSource_, USp_.primitiveFieldRef());

    vsm().mapToSurface(pnSource_, pnSp_.primitiveFieldRef());

    DynamicList<spatialSpreading::sourceInput> mappedImpacts;
    if (spatialSpreadingPtr_ && !impacts_.empty())
    {
        const labelUList& meshFaces = regionMesh().faceLabels();
        Map<label> areaFaceByMeshFace(meshFaces.size());
        forAll(meshFaces, areaFacei)
        {
            areaFaceByMeshFace.insert(meshFaces[areaFacei], areaFacei);
        }

        const polyBoundaryMesh& patches = primaryMesh().boundaryMesh();
        for (const parcelImpact& impact : impacts_)
        {
            if
            (
                impact.patchi < 0
             || impact.patchi >= patches.size()
             || impact.facei < 0
             || impact.facei >= patches[impact.patchi].size()
            )
            {
                FatalErrorInFunction
                    << "Invalid parcel impact patch/face "
                    << impact.patchi << ':' << impact.facei
                    << abort(FatalError);
            }

            const label meshFace =
                patches[impact.patchi].start() + impact.facei;
            if (!areaFaceByMeshFace.contains(meshFace))
            {
                continue;
            }

            const label areaFace = areaFaceByMeshFace.at(meshFace);
            spatialSpreading::sourceInput source;
            source.face = areaFace;
            source.mass = impact.mass;
            source.momentum = impact.momentum;
            source.normalImpulse = impact.normalImpulse;
            source.diameter = impact.diameter;
            mappedImpacts.append(source);
        }
    }

    if (spatialSpreadingPtr_)
    {
        // Map diameter and compute the mass-weighted mean only when required
        auto tDiameterSp = areaScalarField::New
        (
            "diameterSp",
            IOobjectOption::NO_REGISTER,
            regionMesh(),
            dimensionedScalar(dimLength, Zero)
        );
        areaScalarField& diameterSp = tDiameterSp.ref();

        if (spatialSpreadingPtr_->needsDiameter())
        {
            vsm().mapToSurface
            (
                diameterSource_,
                diameterSp.primitiveFieldRef()
            );

            forAll(diameterSp, i)
            {
                const scalar mass = rhoSp_[i];
                if (Foam::mag(mass) > VSMALL)
                {
                    diameterSp[i] /= mass;
                }
                else
                {
                    diameterSp[i] = 0;
                }
            }

            diameterSp.correctBoundaryConditions();
        }

        // Apply identical spatial weights to every source channel.
        if (mappedImpacts.empty())
        {
            spatialSpreadingPtr_->spread
            (
                rhoSp_, USp_, pnSp_, diameterSp
            );
        }
        else
        {
            spatialSpreadingPtr_->spread
            (
                rhoSp_, USp_, pnSp_, diameterSp, mappedImpacts
            );
        }
    }

    auto tFilmMassSp = areaScalarField::New
    (
        "filmMassSp",
        IOobjectOption::NO_REGISTER,
        regionMesh(),
        dimensionedScalar(dimVelocity, Zero)
    );
    areaScalarField& filmMassSp = tFilmMassSp.ref();

    auto tFilmUSp = areaVectorField::New
    (
        "filmUSp",
        IOobjectOption::NO_REGISTER,
        regionMesh(),
        dimensionedVector(sqr(dimVelocity), Zero)
    );
    areaVectorField& filmUSp = tFilmUSp.ref();

    auto tFilmPnSp = areaScalarField::New
    (
        "filmPnSp",
        IOobjectOption::NO_REGISTER,
        regionMesh(),
        dimensionedScalar(dimPressure, Zero)
    );
    areaScalarField& filmPnSp = tFilmPnSp.ref();

    vsm().mapToSurface
    (
        filmMassSource_,
        filmMassSp.primitiveFieldRef()
    );
    vsm().mapToSurface
    (
        filmMomentumSource_,
        filmUSp.primitiveFieldRef()
    );
    vsm().mapToSurface
    (
        filmPnSource_,
        filmPnSp.primitiveFieldRef()
    );

    rhoSp_ += filmMassSp;
    USp_ += filmUSp;
    pnSp_ += filmPnSp;

    // Always force reset the diameter accumulator for the next time step
    diameterSource_ == Zero;

    // Calculate rate per area
    rhoSp_.primitiveFieldRef() *= rAreaDeltaT/rho_;
    USp_.primitiveFieldRef() *= rAreaDeltaT/rho_;
    pnSp_.primitiveFieldRef() *= rAreaDeltaT;

    rhoSp_.relax();
    pnSp_.relax();
    USp_.relax();
}


void liquidFilmModel::postEvolveRegion()
{
    availableMass_ = (h() - h0_)*rho()*regionMesh().S();
    injection_.correct(availableMass_, cloudMassTrans_, cloudDiameterTrans_);
    liquidFilmBase::postEvolveRegion();
    impacts_.clear();
}


void liquidFilmModel::info()
{
    Info<< "\nSurface film: " << type() << " on patch: ";

    for (const label patchi : this->primaryPatchIDs())
    {
        Info<< ' ' << patchi;
    }
    Info<< endl;

    const DimensionedField<scalar, areaMesh>& sf = regionMesh().S();

    Info<< indent << "min/max(mag(Uf))    = "
        << gMinMaxMag(Uf_.field()) << nl
        << indent << "min/max(delta)     = "
        << gMinMax(h_.field()) << nl
        << indent << "coverage           = "
        << gSum(alpha()().field()*mag(sf.field()))/gSumMag(sf.field()) <<  nl
        << indent << "total mass         = "
        << gSum(availableMass_) << nl;

    injection_.info(Info);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
