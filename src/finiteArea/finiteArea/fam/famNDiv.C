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

#include "famNDiv.H"
#include "fvMesh.H"
#include "faMatrices.H"
#include "faConvectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<faMatrix<Type>>
ndiv
(
    const edgeScalarField& flux,
    GeometricField<Type, fvPatchField, areaMesh>& vf,
    const word& name
)
{
    // Implicit normal-divergence.
    //
    // The explicit counterpart fac::ndiv projects the divergence onto the
    // surface normal, result = n*(n & div). That projection is the tensor
    // P = (n (x) n) applied to a Type field. faMatrix<Type> (like
    // fvMatrix<Type>) carries *scalar* coefficients that act identically on
    // every component, so P cannot be represented implicitly without a
    // block-coupled matrix (which faMatrix does not provide): applying n & (.)
    // couples the components. We therefore return the unprojected implicit
    // divergence (identical to upstream finiteArea), leaving the normal
    // constraint to be enforced explicitly / via the field boundary
    // conditions. Do not "project" the scalar coefficients here - that would
    // be numerically incorrect.
    return fa::convectionScheme<Type>::New
    (
        vf.mesh(),
        flux,
        vf.mesh().divScheme(name)
    ).ref().famDiv(flux, vf);
}


template<class Type>
tmp<faMatrix<Type>>
ndiv
(
    const tmp<edgeScalarField>& tflux,
    GeometricField<Type, fvPatchField, areaMesh>& vf,
    const word& name
)
{
    tmp<faMatrix<Type>> Div(fam::ndiv(tflux(), vf, name));
    tflux.clear();

    return Div;
}


template<class Type>
tmp<faMatrix<Type>>
ndiv
(
    const edgeScalarField& flux,
    GeometricField<Type, fvPatchField, areaMesh>& vf
)
{
    return fam::ndiv(flux, vf, "div("+flux.name()+','+vf.name()+')');
}


template<class Type>
tmp<faMatrix<Type>>
ndiv
(
    const tmp<edgeScalarField>& tflux,
    GeometricField<Type, fvPatchField, areaMesh>& vf
)
{
    tmp<faMatrix<Type>> Div(fam::ndiv(tflux(), vf));
    tflux.clear();

    return Div;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
