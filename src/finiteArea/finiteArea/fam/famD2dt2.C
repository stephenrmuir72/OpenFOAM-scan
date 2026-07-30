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

#include "areaFields.H"
#include "faMatrices.H"
#include "fvmD2dt2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// The reduced-order faMesh is an fvMesh (areaMesh == volMesh,
// faMatrix == fvMatrix) and the d2dt2 discretisation is purely temporal, so
// the finite-area second time derivative is identical to the finite-volume
// one. Forward to fvm::d2dt2 rather than duplicate the scheme hierarchy.

template<class Type>
tmp<faMatrix<Type>> d2dt2
(
    const GeometricField<Type, fvPatchField, areaMesh>& vf
)
{
    return fvm::d2dt2(vf);
}


template<class Type>
tmp<faMatrix<Type>> d2dt2
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, areaMesh>& vf
)
{
    return fvm::d2dt2(rho, vf);
}


template<class Type>
tmp<faMatrix<Type>> d2dt2
(
    const areaScalarField& rho,
    const GeometricField<Type, fvPatchField, areaMesh>& vf
)
{
    return fvm::d2dt2(rho, vf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
