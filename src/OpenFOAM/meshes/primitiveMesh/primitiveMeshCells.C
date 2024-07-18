/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2024 OpenCFD Ltd.
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

#include "primitiveMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::primitiveMesh::calcCells
(
    cellList& cellFaceAddr,
    const labelUList& own,
    const labelUList& nei,
    const label inNCells,
    const faceList& faces
)
{
    label nCells = inNCells;

    if (nCells == -1)
    {
        nCells = -1;

        forAll(own, facei)
        {
            nCells = max(nCells, own[facei]);
        }
        nCells++;
    }


    // Determine nTopological from faces
    label nTopologicalD = 3;
    for (const face& f : faces)
    {
        if (f.size() == 1)
        {
            nTopologicalD = 1;
            break;
        }
        else if (f.size() == 2)
        {
            nTopologicalD = 2;
            break;
        }
    }


    // 1. Count number of faces per cell

    labelList ncf(nCells, Zero);

    forAll(own, facei)
    {
        ncf[own[facei]]++;
    }

    forAll(nei, facei)
    {
        if (nei[facei] >= 0)
        {
            ncf[nei[facei]]++;
        }
    }

    // Create the storage
    cellFaceAddr.setSize(ncf.size());


    // 2. Size and fill cellFaceAddr

    forAll(cellFaceAddr, celli)
    {
        cellFaceAddr[celli].setSize(ncf[celli]);
    }
    ncf = 0;

    forAll(own, facei)
    {
        label celli = own[facei];

        cellFaceAddr[celli][ncf[celli]++] = facei;
    }

    forAll(nei, facei)
    {
        label celli = nei[facei];

        if (celli >= 0)
        {
            cellFaceAddr[celli][ncf[celli]++] = facei;
        }
    }


    if (nTopologicalD == 2)
    {
        // Make sure to sort faces in order so we can circulate the outside
        DynamicList<label> oldToNew;
        forAll(cellFaceAddr, celli)
        {
            auto& c = cellFaceAddr[celli];

            oldToNew.clear();
            oldToNew.append(0);
            label pointi = faces[c[oldToNew.last()]][1];

            // Find next face. TBD: could just bubble sort to avoid allocations
            for (label i = 1; i < c.size(); i++)
            {
                label index = -1;
                label nextPointi = -1;

                for (label j = 1; j < c.size(); j++)
                {
                    if (j != oldToNew.last())
                    {
                        const label facei = c[j];
                        const face& f = faces[facei];

                        if (f[0] == pointi)
                        {
                            index = j;
                            nextPointi = f[1];
                        }
                        else if (f[1] == pointi)
                        {
                            index = j;
                            nextPointi = f[0];
                        }
                    }
                }

                if (index == -1)
                {
                    FatalErrorInFunction << "On cell " << celli
                        << " the faces " << flatOutput(c)
                        << " with vertices "
                        << flatOutput(UIndirectList<face>(faces, c))
                        << " do not form a loop" << exit(FatalError);
                }

                // Step to next face and point on face
                oldToNew.append(index);
                pointi = nextPointi;
            }

            // Re-order
            inplaceReorder(oldToNew, c);
        }
    }
}


void Foam::primitiveMesh::calcCells() const
{
    // Loop through faceCells and mark up neighbours

    if (debug)
    {
        Pout<< "primitiveMesh::calcCells() : calculating cells"
            << endl;
    }

    // It is an error to attempt to recalculate cells
    // if the pointer is already set
    if (cfPtr_)
    {
        FatalErrorInFunction
            << "cells already calculated"
            << abort(FatalError);
    }
    else
    {
        // Create the storage
        cfPtr_ = std::make_unique<cellList>(nCells());
        auto& cellFaceAddr = *cfPtr_;

        calcCells
        (
            cellFaceAddr,
            faceOwner(),
            faceNeighbour(),
            nCells(),
            faces()
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::cellList& Foam::primitiveMesh::cells() const
{
    if (!cfPtr_)
    {
        calcCells();
    }

    return *cfPtr_;
}


Foam::face Foam::primitiveMesh::loop(const label celli) const
{
    const auto& fcs = faces();
    const auto& owner = faceOwner();
    const auto& cFaces = cells()[celli];

    for (const label facei : cFaces)
    {
        if (fcs[facei].size() != 2)
        {
            FatalErrorInFunction << "Only applicable for 2D faces"
                << " cell:" << celli
                << " with faces " << cFaces
                << " with vertices " << UIndirectList<face>(fcs, cFaces)
                << exit(FatalError);
        }
    }

    // Storage for new face
    face newF(cFaces.size());
    newF.UList<label>::operator=(-1);
    label fp = 0;

    // Start with zero'th face
    const label facei = cFaces[0];
    edge e(fcs[facei].edge(0));
    if (owner[facei] != celli)
    {
        std::swap(e[0], e[1]);
    }
    newF[fp++] = e[0];
    newF[fp++] = e[1];
    label nextPointi = e[1];
    label prevFacei = facei;

    for (label i = 1; i < cFaces.size(); i++)
    {
        // Find nextPointi in other faces (there should be only one other face)
        for (label j = 1; j < cFaces.size(); j++)
        {
            const label facei = cFaces[j];

            if (facei != prevFacei)
            {
                edge e(fcs[facei].edge(0));
                if (owner[facei] != celli)
                {
                    std::swap(e[0], e[1]);
                }

                if (e[0] == nextPointi)
                {
                    newF[fp++] = e[1];
                    nextPointi = e[1];
                    prevFacei = facei;
                    break;
                }
            }
        }
    }

    if (newF.found(-1))
    {
        FatalErrorInFunction << "Did not find loop on"
            << " cell:" << celli
            << " with faces " << cFaces
            << " with vertices " << UIndirectList<face>(fcs, cFaces)
            << exit(FatalError);
    }

    return newF;
}


// ************************************************************************* //
