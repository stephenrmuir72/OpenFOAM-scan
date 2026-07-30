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

#include "areaFields.H"
#include "faMatrices.H"
#include "fvmDdt.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// The reduced-order faMesh is an fvMesh (areaMesh == volMesh,
// faMatrix == fvMatrix) and the ddt discretisation is purely temporal, so the
// finite-area first time derivative is identical to the finite-volume one.
// Forward to fvm::ddt rather than duplicate the scheme hierarchy.

template<class Type>
tmp<faMatrix<Type>>
ddt
(
    const GeometricField<Type, fvPatchField, areaMesh>& vf
)
{
    return fvm::ddt(vf);
}


template<class Type>
tmp<faMatrix<Type>>
ddt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, areaMesh>& vf
)
{
    return fvm::ddt(rho, vf);
}


template<class Type>
tmp<faMatrix<Type>>
ddt
(
    const areaScalarField& rho,
    const GeometricField<Type, fvPatchField, areaMesh>& vf
)
{
    return fvm::ddt(rho, vf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
