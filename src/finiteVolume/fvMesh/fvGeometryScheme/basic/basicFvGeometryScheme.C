/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "basicFvGeometryScheme.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "primitiveMeshTools.H"
#include "syncTools.H"
#include "triPointRef.H"
#include "wedgePolyPatch.H"
#include "labelIOList.H"
#include "vectorIOField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(basicFvGeometryScheme, 0);
    addToRunTimeSelectionTable(fvGeometryScheme, basicFvGeometryScheme, dict);
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Define an area-weighted normal for three points (defining a triangle)
// (p0, p1, p2) are the base, first, second points respectively
//
// From the original Tukovic code:
//
//   vector n = (d1^d2)/mag(d1^d2);
//   scalar sinAlpha = mag(d1^d2)/(mag(d1)*mag(d2));
//   scalar w = sinAlpha/(mag(d1)*mag(d2));
//
// ie, normal weighted by area, sine angle and inverse distance squared.
// - area : larger weight for larger areas
// - sin  : lower weight for narrow angles (eg, shards)
// - inv distance squared : lower weights for distant points
//
// The above refactored, with 0.5 for area:
//
//   (d1 ^ d2) / (2 * magSqr(d1) * magSqr(d2))

static inline vector areaInvDistSqrWeightedNormal
(
    const vector& a,
    const vector& b
)
{
    // TBD: Use volatile to avoid aggressive branch optimization
    const scalar s(2*magSqr(a)*magSqr(b));

    return s < ROOTVSMALL ? Zero : (a ^ b) / s;
}


// The area normal for the face dual (around base-point)
// described by the right/left edge points and the centre point
//
// The adjustment for 1/2 edge point (the dual point) is done internally
static inline vector areaInvDistSqrWeightedNormalDualEdge
(
    const point& basePoint,
    const point& rightPoint,
    const point& leftPoint,
    const point& centrePoint
)
{
    const vector mid(centrePoint - basePoint);

    return
    (
        areaInvDistSqrWeightedNormal
        (
            0.5*(rightPoint - basePoint),   // vector to right 1/2 edge
            mid
        )
      + areaInvDistSqrWeightedNormal
        (
            mid,
            0.5*(leftPoint - basePoint)     // vector to left 1/2 edge
        )
    );
}


// Weighted area normal for the face triangle (around base-point)
// to the edge points and the centre point.
// Used for example for area-weighted edge normals
// ---------------------------------------------------------------------------
//         own          |
//          *           | Don't bother trying to split the triangle into
//         /.\          | sub-triangles. Planar anyhow and skew weighting
//        / . \         | is less relevant here.
//       /  .  \        |
//      /   .   \       | Note: need negative for neighbour side orientation.
//  e0 *----:----* e1   |
// ---------------------------------------------------------------------------

static inline vector areaInvDistSqrWeightedNormalFaceTriangle
(
    const linePointRef& edgeLine,
    const point& faceCentre
)
{
    return
    (
        areaInvDistSqrWeightedNormal
        (
            (edgeLine.a() - faceCentre),  // From centre to right edge
            (edgeLine.b() - faceCentre)   // From centre to left edge
        )
    );
}


// Simple area normal calculation for specified edge vector and face centre
// The face centre comes last since it is less accurate
static inline vector areaNormalFaceTriangle
(
    const linePointRef& edgeLine,
    const point& faceCentre
)
{
    return triPointRef::areaNormal(edgeLine.a(), edgeLine.b(), faceCentre);
}

} // End namespace Foam


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

//- Calculate the 'Le' vector from faceCentre to edge centre
//  using the edge normal to correct for curvature
//
//  Normalise and rescaled to the edge length
Foam::vector Foam::basicFvGeometryScheme::calcLeVector
(
    const point& faceCentre,
    const linePointRef& edgeLine,
    const vector& edgeNormal   // (unit or area normal)
)
{
    const vector centreToEdge(edgeLine.centre() - faceCentre);

    vector leVector(edgeLine.vec() ^ edgeNormal);

    scalar s(mag(leVector));

    if (s < ROOTVSMALL)
    {
        // The calculated edgeNormal somehow degenerate and thus a
        // bad cross-product?
        // Revert to basic centre -> edge

        leVector = centreToEdge;
        leVector.removeCollinear(edgeLine.unitVec());
        s = mag(leVector);

        if (s < ROOTVSMALL)
        {
            // Unlikely that this should happen
            return Zero;
        }

        leVector *= edgeLine.mag()/s;
    }
    else
    {
        // The additional orientation is probably unnecessary
        leVector *= edgeLine.mag()/s * sign(leVector & centreToEdge);
    }

    return leVector;
}

} // End namespace Foam


namespace Foam
{

// Impose a minimum (edge) vector length
// - disallow any mag(vec) < SMALL
static inline void imposeMinVectorLength(vectorField& fld)
{
    // If it is small, no orientation information reasonably possible
    // so use a vector(1,1,1) * sqrt(1/3)*min-length

    // sqrt(1/3) = 0.5773502691896257, but slightly rounded down
    const vector minVector(vector::uniform(0.57735*SMALL));
    const scalar minLenSqr(SMALL*SMALL);

    for (vector& v : fld)
    {
        if (v.magSqr() < minLenSqr)
        {
            v = minVector;
        }
    }
}

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::basicFvGeometryScheme::calcPointAreaNormals
(
    vectorField& result
) const
{
    // See finiteArea/faMesh/faMeshDemandDrivenData.C

    DebugInFunction
        << "Calculating pointAreaNormals : face dual method" << endl;

    result.resize_nocopy(mesh_.nPoints());
    result = Zero;

    const auto& cells = mesh_.cells();
    const auto& points = mesh_.points();
    const auto& faces = mesh_.faces();
    //const auto& cc = mesh_.cellCentres();

    // Loop over all cells
    //face f;
    forAll(cells, celli)
    {
        // Collect vertices as a face
        const face f(mesh_.loop(celli));

        if (uniqueSort(f).size() != f.size())
        {
            const auto& cFaces = cells[celli];
            FatalErrorInFunction
                << "problem cell:" << celli
                << " faces:" << cFaces
                << " verts:" << UIndirectList<face>(faces, cFaces)
                << " sorted verts:" << f
                << exit(FatalError);
        }

        //Pout<< "cell:" << celli
        //    << " verts:" << flatOutput(f)
        //    << " points:" << flatOutput(UIndirectList<point>(points, f))
        //    << " normal:" << f.unitNormal(points)
        //    << endl;

        //- Use centroid
        //point centrePoint(cc[celli]);
        //- Use average
        point centrePoint(Zero);
        for (const label pointi : f)
        {
            centrePoint += points[pointi];
        }
        centrePoint /= f.size();

        forAll(f, i)
        {
            const label pt0 = f.thisLabel(i);    // base
            result[pt0] +=
                areaInvDistSqrWeightedNormalDualEdge
                (
                    points[pt0],             // base
                    points[f.nextLabel(i)],  // right
                    points[f.prevLabel(i)],  // left
                    centrePoint
                );
        }
    }

    // Parallel: sum the locally-accumulated contributions of each shared point
    // across point coupling.
    //
    // A reduced-order manifold that is a boundary region of a parent volume
    // mesh (the finiteArea faMesh) publishes in its registry the parent-mesh
    // point label of each of its points (a labelList named as below). Its OWN
    // point coupling is degenerate: the edges-as-faces processor patches
    // collapse the global-point addressing, so a shared point is treated as
    // coupled to every processor-boundary point and the summed normal is
    // corrupted at multiply-shared points (e.g. the poles of a decomposed
    // sphere, shared by non-adjacent ranks at a domain corner). When that map
    // and a parent mesh are available, route the sum through the PARENT mesh's
    // (correct) global point coupling, which resolves both adjacent-rank and
    // multiply-shared points in a single collective (the parent-global analogue
    // of the legacy faMesh::calcPointAreaNormals processor point exchange +
    // sharedPoints combineReduce). Otherwise fall back to this mesh's own
    // coupling: correct for a genuine fvMesh, and a no-op in serial.
    //
    // Name kept in sync with faMesh::registerParentPointMap (finiteArea). The
    // seam is a registry lookup (not a virtual/data member on fvMesh) so it does
    // not perturb the fvMesh ABI for the many fvMesh subclasses.
    const labelIOList* parentMapPtr =
        mesh_.thisDb().cfindObject<labelIOList>("pointAreaNormalsParentMap");

    refPtr<fvMesh> parentRef;
    if (parentMapPtr && parentMapPtr->size() == result.size())
    {
        parentRef = mesh_.parentMesh();
    }

    if (parentRef.good())
    {
        syncTools::syncPointList
        (
            parentRef(),
            *parentMapPtr,          // reduced-order point -> parent point
            result,
            plusEqOp<point>(),
            vector::zero,
            mapDistribute::transform()
        );
    }
    else
    {
        syncTools::syncPointList
        (
            mesh_,
            result,
            plusEqOp<point>(),
            vector::zero,
            mapDistribute::transform()
        );
    }

    // TBD: wedge correction?
    {
        for (const auto& fpp : mesh_.boundary())
        {
            if (isA<wedgePolyPatch>(fpp.patch()))
            {
                WarningInFunction<< "Special wedge handling not supported"
                    << endl;
            }
        }
    }

    // TBD: guarantee point normal perpendicular to any coupled fvMesh.
    // Debug-only: this fires on every point-normal (re)build; keep it out of
    // normal run output until the coupled-perpendicular correction is added.
    if (debug)
    {
        WarningInFunction
            << "PointNormals not perpendicular to any coupled fvPatch."
            << endl;
    }

//for (label facei = mesh_.nInternalFaces(); facei < mesh_.nFaces(); facei++)
//{
//    const face& f = faces[facei];
//    forAll(f, i)
//    {
//        result[f[i]].removeCollinear(vector(0, 0, 1));
//    }
//}

    // Boundary (contact-line) point-normal correction seam.
    //
    // A reduced-order manifold coupled to a parent volume mesh (the finiteArea
    // faMesh) may require its boundary point normals to be made perpendicular to
    // a neighbouring wall polyPatch, so that the finite-area curvature -- the
    // surface-tension restoring term in free-surface tracking -- is correct at
    // the contact line. That coupling (neighbour polyPatch face normals, and
    // which patches opt in) is faMesh-specific and cannot be derived by this
    // generic scheme, so faMesh publishes, per reduced-order point, the
    // accumulated neighbour-face normal whose component is to be removed (Zero
    // elsewhere), refreshed before each geometry (re)build.
    //
    // It is applied HERE -- inside the point-normal build, hence before the
    // eager edge derivation in faGeometryScheme::movePoints (edgeAreaNormals,
    // Le, rotated edge Sf) and therefore before faceCurvatures -- exactly as the
    // legacy faMesh::calcPointAreaNormals removeCollinear(haloFaceNormals)
    // correction. A registry entry (not an fvMesh member/virtual) keeps the
    // fvMesh ABI untouched, mirroring the "pointAreaNormalsParentMap" seam above.
    //
    // Name kept in sync with faMesh::updatePointNormalsBoundaryCorrection.
    const vectorIOField* pnCorrPtr =
        mesh_.thisDb().cfindObject<vectorIOField>
        (
            "pointAreaNormalsBoundaryCorrection"
        );

    if (pnCorrPtr && pnCorrPtr->size() == result.size())
    {
        const vectorField& pnCorr = *pnCorrPtr;
        forAll(result, pointi)
        {
            if (magSqr(pnCorr[pointi]) > ROOTVSMALL)
            {
                result[pointi].removeCollinear(normalised(pnCorr[pointi]));
            }
        }
    }

    result.normalise();
}


void Foam::basicFvGeometryScheme::calcFaceAreaNormalsInvDistSqrWeighted
(
    vectorField& fDir1
) const
{
    // See finiteArea/faMesh/faMeshDemandDrivenData.C

    DebugInFunction
        << "Calculating faceAreaNormals : face triangle-normal method" << endl;

    const labelUList& owner = mesh_.faceOwner();
    const labelUList& neighbour = mesh_.neighbour();
    const auto& faces = mesh_.faces();
    const auto& pts = mesh_.points();
    const auto& cc = mesh_.cellCentres();

    forAll(owner, facei)
    {
        const face& f = faces[facei];
        const linePointRef edgeLine(pts[f[0]], pts[f[1]]);

        fDir1[facei] = areaInvDistSqrWeightedNormalFaceTriangle
        (
            edgeLine,
            cc[owner[facei]]
        );
    }

    vectorField nbrDir1
    (
        SubList<vector>
        (
            fDir1,
            mesh_.nBoundaryFaces(),
            mesh_.nInternalFaces()
        )
    );
    syncTools::swapBoundaryFaceList(mesh_, nbrDir1);

    // Finish internal faces
    forAll(neighbour, facei)
    {
        const face& f = faces[facei];
        const linePointRef edgeLine(pts[f[1]], pts[f[0]]);

        fDir1[facei] += areaInvDistSqrWeightedNormalFaceTriangle
        (
            edgeLine,
            cc[neighbour[facei]]
        );
    }

    // Finish coupled faces
    for (auto& fpp : mesh_.boundary())
    {
        if (isA<wedgePolyPatch>(fpp.patch()))
        {
            WarningInFunction<< "Special wedge handling not supported"
                << endl;
        }

        const auto* ppPtr = isA<coupledPolyPatch>(fpp.patch());
        if (ppPtr)
        {
            const auto& pp = *ppPtr;
            forAll(pp, i)
            {
                const label bFacei = pp.offset()+i;
                fDir1[pp.start()+i] += nbrDir1[bFacei];
            }
        }
    }


    // TBD: guarantee point normal perpendicular to any coupled fvMesh.
    // Debug-only: this fires on every point-normal (re)build; keep it out of
    // normal run output until the coupled-perpendicular correction is added.
    if (debug)
    {
        WarningInFunction
            << "PointNormals not perpendicular to any coupled fvPatch."
            << endl;
    }

//for (label facei = mesh_.nInternalFaces(); facei < mesh_.nFaces(); facei++)
//{
//    fDir1[facei].removeCollinear(vector(0, 0, 1));
//}

    fDir1.normalise();
    imposeMinVectorLength(fDir1);
}


void Foam::basicFvGeometryScheme::calcFaceAreaNormalsAvgPoint
(
    const vectorField& pointNormals,
    vectorField& fDir1
) const
{
    DebugInFunction
        << "Calculating faceAreaNormals : average point-normal method" << endl;

    const auto& faces = mesh_.faces();
    const auto& pts = mesh_.points();

    // All faces
    forAll(faces, facei)
    {
        const face& f = faces[facei];
        const linePointRef edgeLine(pts[f[0]], pts[f[1]]);

        // Average of both ends
        fDir1[facei] = (pointNormals[f[0]] + pointNormals[f[1]]);

        fDir1[facei].removeCollinear(edgeLine.unitVec());
        fDir1[facei].normalise();
    }

    imposeMinVectorLength(fDir1);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicFvGeometryScheme::basicFvGeometryScheme
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvGeometryScheme(mesh, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::basicFvGeometryScheme::movePoints()
{
    fvGeometryScheme::movePoints();

    if (debug)
    {
        Pout<< "basicFvGeometryScheme::movePoints() : "
            << "recalculating primitiveMesh centres" << endl;
    }

    // Use lower level to calculate the geometry.
    // Reduced-order (2D) manifold geometry is handled by the dedicated
    // faGeometryScheme (selected automatically for 2D meshes).
    const_cast<fvMesh&>(mesh_).primitiveMesh::updateGeom();
}


Foam::tmp<Foam::surfaceScalarField> Foam::basicFvGeometryScheme::weights() const
{
    if (debug)
    {
        Pout<< "basicFvGeometryScheme::weights() : "
            << "Constructing weighting factors for face interpolation"
            << endl;
    }

    auto tweights =
        tmp<surfaceScalarField>::New
        (
            IOobject
            (
                "weights",
                mesh_.pointsInstance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh_,
            dimless
        );

    auto& weights = tweights.ref();
    weights.setOriented();

    // Set local references to mesh data
    // Note that we should not use fvMesh sliced fields at this point yet
    // since this causes a loop when generating weighting factors in
    // coupledFvPatchField evaluation phase
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();

    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& C = mesh_.cellCentres();
    const vectorField& Sf = mesh_.faceAreas();

    // ... and reference to the internal field of the weighting factors
    scalarField& w = weights.primitiveFieldRef();

    forAll(owner, facei)
    {
        // Note: mag in the dot-product.
        // For all valid meshes, the non-orthogonality will be less than
        // 90 deg and the dot-product will be positive.  For invalid
        // meshes (d & s <= 0), this will stabilise the calculation
        // but the result will be poor.
        scalar SfdOwn = mag(Sf[facei] & (Cf[facei] - C[owner[facei]]));
        scalar SfdNei = mag(Sf[facei] & (C[neighbour[facei]] - Cf[facei]));

        if (mag(SfdOwn + SfdNei) > ROOTVSMALL)
        {
            w[facei] = SfdNei/(SfdOwn + SfdNei);
        }
        else
        {
            w[facei] = 0.5;
        }
    }

    auto& wBf = weights.boundaryFieldRef();

    forAll(mesh_.boundary(), patchi)
    {
        mesh_.boundary()[patchi].makeWeights(wBf[patchi]);
    }

    if (debug)
    {
        Pout<< "basicFvGeometryScheme::weights : "
            << "Finished constructing weighting factors for face interpolation"
            << endl;
    }
    return tweights;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::basicFvGeometryScheme::deltaCoeffs() const
{
    if (debug)
    {
        Pout<< "basicFvGeometryScheme::deltaCoeffs() : "
            << "Constructing differencing factors array for face gradient"
            << endl;
    }

    // Force the construction of the weighting factors
    // needed to make sure deltaCoeffs are calculated for parallel runs.
    (void)mesh_.weights();

    auto tdeltaCoeffs =
        tmp<surfaceScalarField>::New
        (
            IOobject
            (
                "deltaCoeffs",
                mesh_.pointsInstance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh_,
            dimless/dimLength
        );

    auto& deltaCoeffs = tdeltaCoeffs.ref();
    deltaCoeffs.setOriented();


    // Set local references to mesh data
    const vectorField& C = mesh_.cellCentres();
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();

    forAll(owner, facei)
    {
        deltaCoeffs[facei] = 1.0/mag(C[neighbour[facei]] - C[owner[facei]]);
    }

    auto& deltaCoeffsBf = deltaCoeffs.boundaryFieldRef();

    forAll(deltaCoeffsBf, patchi)
    {
        const fvPatch& p = mesh_.boundary()[patchi];
        deltaCoeffsBf[patchi] = 1.0/mag(p.delta());

        // Optionally correct
        p.makeDeltaCoeffs(deltaCoeffsBf[patchi]);
    }

    return tdeltaCoeffs;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::basicFvGeometryScheme::nonOrthDeltaCoeffs() const
{
    if (debug)
    {
        Pout<< "basicFvGeometryScheme::nonOrthDeltaCoeffs() : "
            << "Constructing differencing factors array for face gradient"
            << endl;
    }

    // Force the construction of the weighting factors
    // needed to make sure deltaCoeffs are calculated for parallel runs.
    weights();

    auto tnonOrthDeltaCoeffs =
        tmp<surfaceScalarField>::New
        (
            IOobject
            (
                "nonOrthDeltaCoeffs",
                mesh_.pointsInstance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh_,
            dimless/dimLength
        );

    auto& nonOrthDeltaCoeffs = tnonOrthDeltaCoeffs.ref();
    nonOrthDeltaCoeffs.setOriented();


    // Set local references to mesh data
    const volVectorField& C = mesh_.C();
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();

    forAll(owner, facei)
    {
        vector delta = C[neighbour[facei]] - C[owner[facei]];
        vector unitArea = Sf[facei]/magSf[facei];

        // Standard cell-centre distance form
        //NonOrthDeltaCoeffs[facei] = (unitArea & delta)/magSqr(delta);

        // Slightly under-relaxed form
        //NonOrthDeltaCoeffs[facei] = 1.0/mag(delta);

        // More under-relaxed form
        //NonOrthDeltaCoeffs[facei] = 1.0/(mag(unitArea & delta) + VSMALL);

        // Stabilised form for bad meshes
        nonOrthDeltaCoeffs[facei] = 1.0/max(unitArea & delta, 0.05*mag(delta));
    }

    auto& nonOrthDeltaCoeffsBf = nonOrthDeltaCoeffs.boundaryFieldRef();

    forAll(nonOrthDeltaCoeffsBf, patchi)
    {
        fvsPatchScalarField& patchDeltaCoeffs = nonOrthDeltaCoeffsBf[patchi];

        const fvPatch& p = patchDeltaCoeffs.patch();

        const vectorField patchDeltas(mesh_.boundary()[patchi].delta());

        forAll(p, patchFacei)
        {
            vector unitArea =
                Sf.boundaryField()[patchi][patchFacei]
               /magSf.boundaryField()[patchi][patchFacei];

            const vector& delta = patchDeltas[patchFacei];

            patchDeltaCoeffs[patchFacei] =
                1.0/max(unitArea & delta, 0.05*mag(delta));
        }

        // Optionally correct
        p.makeNonOrthoDeltaCoeffs(patchDeltaCoeffs);
    }
    return tnonOrthDeltaCoeffs;
}


Foam::tmp<Foam::surfaceVectorField>
Foam::basicFvGeometryScheme::nonOrthCorrectionVectors() const
{
    if (debug)
    {
        Pout<< "surfaceInterpolation::makeNonOrthCorrectionVectors() : "
            << "Constructing non-orthogonal correction vectors"
            << endl;
    }

    auto tnonOrthCorrectionVectors =
        tmp<surfaceVectorField>::New
        (
            IOobject
            (
                "nonOrthCorrectionVectors",
                mesh_.pointsInstance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh_,
            dimless
        );

    auto& corrVecs = tnonOrthCorrectionVectors.ref();
    corrVecs.setOriented();

    // Set local references to mesh data
    const volVectorField& C = mesh_.C();
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();
    tmp<surfaceScalarField> tNonOrthDeltaCoeffs(nonOrthDeltaCoeffs());
    const surfaceScalarField& NonOrthDeltaCoeffs = tNonOrthDeltaCoeffs();

    forAll(owner, facei)
    {
        vector unitArea(Sf[facei]/magSf[facei]);
        vector delta(C[neighbour[facei]] - C[owner[facei]]);

        corrVecs[facei] = unitArea - delta*NonOrthDeltaCoeffs[facei];
    }

    // Boundary correction vectors set to zero for boundary patches
    // and calculated consistently with internal corrections for
    // coupled patches

    auto& corrVecsBf = corrVecs.boundaryFieldRef();

    forAll(corrVecsBf, patchi)
    {
        fvsPatchVectorField& patchCorrVecs = corrVecsBf[patchi];

        const fvPatch& p = patchCorrVecs.patch();

        if (!patchCorrVecs.coupled())
        {
            patchCorrVecs = Zero;
        }
        else
        {
            const auto& patchNonOrthDeltaCoeffs =
                NonOrthDeltaCoeffs.boundaryField()[patchi];

            const vectorField patchDeltas(mesh_.boundary()[patchi].delta());

            forAll(p, patchFacei)
            {
                vector unitArea =
                    Sf.boundaryField()[patchi][patchFacei]
                   /magSf.boundaryField()[patchi][patchFacei];

                const vector& delta = patchDeltas[patchFacei];

                patchCorrVecs[patchFacei] =
                    unitArea - delta*patchNonOrthDeltaCoeffs[patchFacei];
            }
        }

        // Optionally correct
        p.makeNonOrthoCorrVectors(patchCorrVecs);
    }

    if (debug)
    {
        Pout<< "surfaceInterpolation::makeNonOrthCorrectionVectors() : "
            << "Finished constructing non-orthogonal correction vectors"
            << endl;
    }
    return tnonOrthCorrectionVectors;
}


bool Foam::basicFvGeometryScheme::updateGeom
(
    const pointField& points,
    const refPtr<pointField>& oldPoints,
    pointField& faceCentres,
    vectorField& faceAreas,
    pointField& cellCentres,
    scalarField& cellVolumes
) const
{
    primitiveMeshTools::makeFaceCentresAndAreas
    (
        mesh_,
        points,
        faceCentres,
        faceAreas
    );

    // Cell direction (element normal) is only meaningful for reduced-order
    // (2D) meshes; leave empty otherwise so makeCellCentresAndVols() skips it.
    vectorField cellDir0;
    if (mesh_.nTopologicalD() == 2)
    {
        cellDir0.setSize(mesh_.nCells(), Zero);
    }

    primitiveMeshTools::makeCellCentresAndVols
    (
        mesh_,
        faceCentres,
        faceAreas,
        cellCentres,
        cellVolumes,
        cellDir0
    );

    // Assume something has changed
    return true;
}


// ************************************************************************* //
