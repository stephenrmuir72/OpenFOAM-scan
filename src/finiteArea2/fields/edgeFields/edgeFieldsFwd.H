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

InNamespace
    Foam

Description
    Forwards for edge field types

\*---------------------------------------------------------------------------*/

#ifndef Foam_edgeFieldsFwd_H
#define Foam_edgeFieldsFwd_H

#include "fieldTypes.H"
//#include "GeometricFieldFwd.H"
#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef surfaceMesh edgeMesh;

// Forward Declarations

//- An edge patch field for a given type
template<class Type>
using faePatchField = fvsPatchField<Type>;

//
//
///// Future? (TBD)
///// //- An (area) edge field for a given type
///// template<class Type>
///// using EdgeField = GeometricField<Type, faePatchField, edgeMesh>;
//
//
//// Typedefs
//
typedef
    GeometricField<scalar, fvsPatchField, edgeMesh>
    edgeScalarField;

typedef
    GeometricField<vector, fvsPatchField, edgeMesh>
    edgeVectorField;

typedef
    GeometricField<sphericalTensor, fvsPatchField, edgeMesh>
    edgeSphericalTensorField;

typedef
    GeometricField<symmTensor, fvsPatchField, edgeMesh>
    edgeSymmTensorField;

typedef
    GeometricField<tensor, fvsPatchField, edgeMesh>
    edgeTensorField;

typedef faePatchField<scalar> faePatchScalarField;
typedef faePatchField<vector> faePatchVectorField;
typedef faePatchField<sphericalTensor> faePatchSphericalTensorField;
typedef faePatchField<symmTensor> faePatchSymmTensorField;
typedef faePatchField<tensor> faePatchTensorField;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
