/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021-2024 OpenCFD Ltd.
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

Description
    Efficient cell-centre calculation using face-addressing, face-centres and
    face-areas.

\*---------------------------------------------------------------------------*/

#include "primitiveMesh.H"
#include "primitiveMeshTools.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::primitiveMesh::calcCellCentresAndVols() const
{
    if (debug)
    {
        Pout<< "primitiveMesh::calcCellCentresAndVols() : "
            << "Calculating cell centres and volumes"
            << endl;
    }

    // These are always calculated in tandem, but only once.
    if (cellCentresPtr_ || cellVolumesPtr_)
    {
        FatalErrorInFunction
            << "Cell centres or volumes already calculated"
            << abort(FatalError);
    }

    // set the accumulated cell centre to zero vector
    cellCentresPtr_ = std::make_unique<vectorField>(nCells(), Zero);
    auto& cellCtrs = *cellCentresPtr_;

    // Initialise cell volumes to 0
    cellVolumesPtr_ = std::make_unique<scalarField>(nCells(), Zero);
    auto& cellVols = *cellVolumesPtr_;

    cellDir0Ptr_ = new vectorField();
    vectorField& cellDir0 = *cellDir0Ptr_;
    if (nTopologicalD() == 2)
    {
        cellDir0.setSize(nCells(), Zero);
    }

    // Make centres and volumes
    primitiveMeshTools::makeCellCentresAndVols
    (
        *this,
        faceCentres(),
        faceAreas(),
        cellCtrs,
        cellVols,
        cellDir0
    );

    if (debug)
    {
        Pout<< "primitiveMesh::calcCellCentresAndVols() : "
            << "Finished calculating cell centres and volumes"
            << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::vectorField& Foam::primitiveMesh::cellCentres() const
{
    if (!cellCentresPtr_)
    {
        //calcCellCentresAndVols();
        const_cast<primitiveMesh&>(*this).updateGeom();
    }

    return *cellCentresPtr_;
}


const Foam::scalarField& Foam::primitiveMesh::cellVolumes() const
{
    if (!cellVolumesPtr_)
    {
        //calcCellCentresAndVols();
        const_cast<primitiveMesh&>(*this).updateGeom();
    }

    return *cellVolumesPtr_;
}


const Foam::vectorField& Foam::primitiveMesh::cellDir0() const
{
    if (!cellDir0Ptr_)
    {
        //calcCellCentresAndVols();
        const_cast<primitiveMesh&>(*this).updateGeom();
    }

    return *cellDir0Ptr_;
}


Foam::tmp<Foam::vectorField> Foam::primitiveMesh::cellAreas() const
{
    return cellVolumes()*cellDir0();
}


const Foam::vectorField& Foam::primitiveMesh::cellDir1() const
{
    if (!cellDir1Ptr_)
    {
        // Externally calculated
        FatalErrorInFunction<< "No calculated tangential vector"
            << exit(FatalError);
    }

    return *cellDir1Ptr_;
}


const Foam::vectorField& Foam::primitiveMesh::pointDir0() const
{
    if (!pointDir0Ptr_)
    {
        // Externally calculated
        FatalErrorInFunction<< "No calculated point normal direction"
            << exit(FatalError);
    }

    return *pointDir0Ptr_;
}


// ************************************************************************* //
