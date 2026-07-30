/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024-2026 OpenCFD Ltd.
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

#include "faGeometryScheme.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "syncTools.H"
#include "coupledPolyPatch.H"
#include "linePointRef.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faGeometryScheme, 0);
    addToRunTimeSelectionTable(fvGeometryScheme, faGeometryScheme, dict);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faGeometryScheme::faGeometryScheme
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    basicFvGeometryScheme(mesh, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faGeometryScheme::movePoints()
{
    // Standard (3D-style) geometry first.
    // This also runs primitiveMesh::updateGeom() via the base class.
    basicFvGeometryScheme::movePoints();

    if (mesh_.nTopologicalD() == 2)
    {
        // In 2D above primitiveMesh will now have the edge-vectors stored
        // in the faceAreas.
        // Correct:
        // - calculate cell direction1 (= normal, since is now a face)
        // - option1: average cell direction1 onto face
        // - option2: average triangle direction1 onto face
        // - re-orient face area to be perpendicular to direction1

        const labelUList& owner = mesh_.faceOwner();
        const auto& faces = mesh_.faces();
        const auto& pts = mesh_.points();
        const auto& cc = mesh_.cellCentres();

        // Calculate point normals from cell normals (cellDir0)
        autoPtr<vectorField> tpointDir0
        (
            autoPtr<vectorField>::New(mesh_.nPoints())
        );
        calcPointAreaNormals(tpointDir0.ref());


        // Average of element normals onto face
        autoPtr<vectorField> tfDir1(autoPtr<vectorField>::New(mesh_.nFaces()));
        auto& fDir1 = tfDir1.ref();


        // Face-tangential normal fDir1 from an average of the point-area
        // normals (calcFaceAreaNormalsAvgPoint). The alternative triangle
        // inverse-distance-weighted construction (calcFaceAreaNormalsInvDist-
        // SqrWeighted) and the average-cell-normal construction are retained as
        // selectable helpers but are not the default; switch here if a future
        // scheme needs them.
        calcFaceAreaNormalsAvgPoint(tpointDir0(), fDir1);


        vectorField& faceAreas =
            const_cast<vectorField&>(mesh_.faceAreas());

        // Restore the raw 2D face areas (the edge vectors p[f[1]]-p[f[0]], see
        // primitiveMeshTools::makeFaceCentresAndAreas) before applying the
        // in-plane rotation below. primitiveMesh::updateGeom() only (re)computes
        // faceAreas when its cache is empty, so a second movePoints() with no
        // intervening clearGeom would otherwise feed the already-rotated areas
        // back into the rotation and corrupt them (non-idempotent): the edge
        // "normal" would be spun out of the manifold tangent plane, collapsing
        // deltaCoeffs/nf on both internal and boundary faces. Re-deriving the
        // base edge vectors here makes the correction idempotent.
        forAll(faceAreas, facei)
        {
            const face& f = faces[facei];
            faceAreas[facei] = pts[f[1]] - pts[f[0]];
        }

        autoPtr<vectorField> tLe(autoPtr<vectorField>::New(faceAreas));
        // Re-calculate edge vectors from the edgeAreaNormals
        {
            vectorField& Le = tLe.ref();

            forAll(owner, facei)
            {
                const face& f = faces[facei];
                const linePointRef edgeLine(pts[f[0]], pts[f[1]]);

                Le[facei] = calcLeVector
                (
                    cc[owner[facei]],   // only used for sign
                    edgeLine,
                    fDir1[facei]        // face-tangential direction
                );
            }
        }


        // Re-do face-area to be perpendicular to the tangential
        // face-directions
        const scalarField magS(mag(faceAreas));

        faceAreas = (faceAreas ^ fDir1);
        faceAreas.normalise();
        faceAreas *= magS;

        // Store on the primitiveMesh level; fvMesh will demand-driven slice it
        // into vol- and surface- fields
        const_cast<fvMesh&>(mesh_).primitiveMesh::resetPrimitives
        (
            autoPtr<vectorField>(),
            autoPtr<vectorField>(),
            std::move(tfDir1),          // face normal from average cell normal
            std::move(tLe),             // edge vector (unnormalised)
            std::move(tpointDir0)       // point normal
        );
    }
}


// ************************************************************************* //
