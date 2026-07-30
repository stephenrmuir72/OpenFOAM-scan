/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2022 OpenCFD Ltd.
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
#include "edgeFields.H"
#include "facGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

inline tmp<areaScalarField> limiter(const areaScalarField& phi)
{
    return phi;
}

inline tmp<areaScalarField> limiter(const areaVectorField& phi)
{
    return magSqr(phi);
}

inline tmp<areaScalarField> limiter(const areaTensorField& phi)
{
    return magSqr(phi);
}

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class NVDweight>
Foam::tmp<Foam::edgeScalarField> Foam::faNVDscheme<Type,NVDweight>::weights
(
    const GeometricField<Type, fvPatchField, areaMesh>& phi
) const
{
    const fvMesh& mesh = this->mesh();

    tmp<edgeScalarField> tWeightingFactors
    (
        new edgeScalarField(mesh.surfaceInterpolation::weights())
    );
    edgeScalarField& weightingFactors = tWeightingFactors.ref();

    scalarField& weights = weightingFactors.primitiveFieldRef();

    tmp<areaScalarField> tvf = limiter(phi);
    const areaScalarField& vf = tvf();

    const areaVectorField gradc(fac::grad(vf));

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const vectorField& n = mesh.faceAreaNormals().internalField();
    const vectorField& c = mesh.areaCentres().internalField();
    const vectorField& Cf = mesh.Cf().internalField();

    forAll(weights, edge)
    {
        vector d(c[neighbour[edge]] - c[owner[edge]]);

        if (edgeFlux_[edge] > 0)
        {
            d.removeCollinear(n[owner[edge]]);
        }
        else
        {
            d.removeCollinear(n[neighbour[edge]]);
        }
        d.normalise();

        // Geodesic P->N distance (owner -> edge centre -> neighbour),
        // formerly edgeInterpolation::lPN(), scales the direction vector.
        d *=
            mag(Cf[edge] - c[owner[edge]])
          + mag(c[neighbour[edge]] - Cf[edge]);

        // No zero-length vectors (after normalise)
        if (d.magSqr() < ROOTSMALL)
        {
            d = vector::uniform(SMALL);
            // OR: d = vector::uniform(0.57735*SMALL);
        }

        weights[edge] =
            this->weight
            (
                weights[edge],
                edgeFlux_[edge],
                vf[owner[edge]],
                vf[neighbour[edge]],
                gradc[owner[edge]],
                gradc[neighbour[edge]],
                d
            );
    }


    auto& bWeights = weightingFactors.boundaryFieldRef();

    forAll(bWeights, patchI)
    {
        if (bWeights[patchI].coupled())
        {
            scalarField& pWeights = bWeights[patchI];

            const scalarField& pEdgeFlux = edgeFlux_.boundaryField()[patchI];

            scalarField pVfP(vf.boundaryField()[patchI].patchInternalField());

            scalarField pVfN(vf.boundaryField()[patchI].patchNeighbourField());

            vectorField pGradcP
            (
                gradc.boundaryField()[patchI].patchInternalField()
            );

            vectorField pGradcN
            (
                gradc.boundaryField()[patchI].patchNeighbourField()
            );

            vectorField CP
            (
                mesh.areaCentres().boundaryField()[patchI].patchInternalField()
            );

            vectorField CN
            (
                mesh.areaCentres().boundaryField()[patchI]
                .patchNeighbourField()
            );

            vectorField nP
            (
                mesh.faceAreaNormals().boundaryField()[patchI]
               .patchInternalField()
            );

            vectorField nN
            (
                mesh.faceAreaNormals().boundaryField()[patchI]
               .patchNeighbourField()
            );

            const vectorField& Cfb = mesh.Cf().boundaryField()[patchI];

            forAll(pWeights, edgeI)
            {
                vector d(CN[edgeI] - CP[edgeI]);

                if (pEdgeFlux[edgeI] > 0)
                {
                    d.removeCollinear(nP[edgeI]);
                }
                else
                {
                    d.removeCollinear(nN[edgeI]);
                }
                d.normalise();

                // Geodesic P->N distance (owner -> edge centre ->
                // neighbour), formerly edgeInterpolation::lPN().
                d *=
                    mag(Cfb[edgeI] - CP[edgeI])
                  + mag(CN[edgeI] - Cfb[edgeI]);

                // No zero-length vectors (after normalise)
                if (d.magSqr() < ROOTSMALL)
                {
                    d = vector::uniform(SMALL);
                    // OR: d = vector::uniform(0.57735*SMALL);
                }

                pWeights[edgeI] =
                    this->weight
                    (
                        pWeights[edgeI],
                        pEdgeFlux[edgeI],
                        pVfP[edgeI],
                        pVfN[edgeI],
                        pGradcP[edgeI],
                        pGradcN[edgeI],
                        d
                    );
            }
        }
    }

    return tWeightingFactors;
}


// ************************************************************************* //
