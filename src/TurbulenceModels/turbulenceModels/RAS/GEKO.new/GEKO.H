/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    GEKO (Generalized k-omega) turbulence model wrapper.

    Mirrors the kOmegaSST wrapper pattern:
    - Inherits from GEKOBase<eddyViscosity<RASModel<BasicTurbulenceModel>>>
    - Calls BasicTurbulenceModel::correctNut() after updating nut_

\*---------------------------------------------------------------------------*/

#ifndef GEKO_H
#define GEKO_H

#include "GEKOBase.H"
#include "RASModel.H"
#include "eddyViscosity.H"

namespace Foam
{
namespace RASModels
{

template<class BasicTurbulenceModel>
class GEKO
:
    public GEKOBase<eddyViscosity<RASModel<BasicTurbulenceModel>>>
{
    typedef GEKOBase<eddyViscosity<RASModel<BasicTurbulenceModel>>> BaseClass;

protected:

    //- Update eddy viscosity and turbulence thermal diffusivity
    virtual void correctNut(const volScalarField& S)
    {
        BaseClass::correctNut(S);
        BasicTurbulenceModel::correctNut();
    }

    //- Update eddy viscosity and turbulence thermal diffusivity
    virtual void correctNut()
    {
        BaseClass::correctNut();
        BasicTurbulenceModel::correctNut();
    }

public:

    TypeName("GEKO");

    GEKO
    (
        const typename BasicTurbulenceModel::alphaField& alpha,
        const typename BasicTurbulenceModel::rhoField& rho,
        const volVectorField& U,
        const surfaceScalarField& alphaRhoPhi,
        const surfaceScalarField& phi,
        const typename BasicTurbulenceModel::transportModel& transport,
        const word& propertiesName = turbulenceModel::propertiesName
    );

    virtual ~GEKO() = default;
};

} // End namespace RASModels
} // End namespace Foam

#endif

// ************************************************************************* //
