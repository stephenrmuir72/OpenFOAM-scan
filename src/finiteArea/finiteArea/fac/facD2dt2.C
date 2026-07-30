/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 Volkswagen AG
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

#include "facD2dt2.H"
#include "fvMesh.H"
#include "fvcD2dt2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fac
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// The reduced-order faMesh is an fvMesh (areaMesh == volMesh) and the d2dt2
// discretisation is purely temporal, so the finite-area explicit second time
// derivative is identical to the finite-volume one. Forward to fvc::d2dt2
// rather than duplicate the scheme hierarchy.
//
// Note: the finite-volume fvc::d2dt2 provides only the field and
// field-scalar-field forms; the (dimensioned<Type>, mesh) and
// (dimensionedScalar, field) forms carried by the legacy finite-area scheme
// had no callers and no fvc counterpart, so they are not reprovided here.

template<class Type>
tmp<GeometricField<Type, fvPatchField, areaMesh>> d2dt2
(
    const GeometricField<Type, fvPatchField, areaMesh>& vf
)
{
    return fvc::d2dt2(vf);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, areaMesh>> d2dt2
(
    const areaScalarField& rho,
    const GeometricField<Type, fvPatchField, areaMesh>& vf
)
{
    return fvc::d2dt2(rho, vf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fac

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
