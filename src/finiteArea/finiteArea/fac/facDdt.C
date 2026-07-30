/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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

#include "facDdt.H"
#include "fvMesh.H"
#include "fvcDdt.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fac
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// The reduced-order faMesh is an fvMesh (areaMesh == volMesh) and the ddt
// discretisation is purely temporal, so the finite-area explicit first time
// derivative is identical to the finite-volume one. Forward to fvc::ddt
// rather than duplicate the scheme hierarchy.

template<class Type>
tmp<GeometricField<Type, fvPatchField, areaMesh>>
ddt
(
    const dimensioned<Type> dt,
    const fvMesh& mesh
)
{
    return fvc::ddt(dt, mesh);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, areaMesh>>
ddt
(
    const GeometricField<Type, fvPatchField, areaMesh>& vf
)
{
    return fvc::ddt(vf);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, areaMesh>>
ddt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, areaMesh>& vf
)
{
    return fvc::ddt(rho, vf);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, areaMesh>>
ddt
(
    const areaScalarField& rho,
    const GeometricField<Type, fvPatchField, areaMesh>& vf
)
{
    return fvc::ddt(rho, vf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fac

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
