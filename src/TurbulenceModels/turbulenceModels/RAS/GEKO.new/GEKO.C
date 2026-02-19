/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
\*---------------------------------------------------------------------------*/

#include "GEKO.H"

namespace Foam
{
namespace RASModels
{

template<class BasicTurbulenceModel>
GEKO<BasicTurbulenceModel>::GEKO
(
    const typename BasicTurbulenceModel::alphaField& alpha,
    const typename BasicTurbulenceModel::rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const typename BasicTurbulenceModel::transportModel& transport,
    const word& propertiesName
)
:
    BaseClass
    (
        "GEKO",
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    )
{}

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
