/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

#include "edgeInterpolation.H"
#include "faMesh.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(edgeInterpolation, 0);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::edgeInterpolation::edgeInterpolation(const faMesh& fam)
:
    faMesh_(fam),
    orthogonal_(false),
    skew_(true)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::surfaceScalarField& Foam::edgeInterpolation::lPN() const
{
    if (!lPNptr_)
    {
        makeLPN();
    }

    return (*lPNptr_);
}


const Foam::surfaceScalarField& Foam::edgeInterpolation::weights() const
{
    if (!weightingFactorsPtr_)
    {
        makeWeights();
    }

    return (*weightingFactorsPtr_);
}


const Foam::surfaceScalarField& Foam::edgeInterpolation::deltaCoeffs() const
{
    if (!differenceFactorsPtr_)
    {
        makeDeltaCoeffs();
    }

    return (*differenceFactorsPtr_);
}


bool Foam::edgeInterpolation::orthogonal() const
{
    if (orthogonal_ == false && !correctionVectorsPtr_)
    {
        makeCorrectionVectors();
    }

    return orthogonal_;
}


const Foam::surfaceVectorField& Foam::edgeInterpolation::correctionVectors() const
{
    if (orthogonal())
    {
        FatalErrorInFunction
            << "cannot return correctionVectors; mesh is orthogonal"
            << abort(FatalError);
    }

    return (*correctionVectorsPtr_);
}


bool Foam::edgeInterpolation::skew() const
{
    if (skew_ == true && !skewCorrectionVectorsPtr_)
    {
        makeSkewCorrectionVectors();
    }

    return skew_;
}


const Foam::surfaceVectorField&
Foam::edgeInterpolation::skewCorrectionVectors() const
{
    if (!skew())
    {
        FatalErrorInFunction
            << "cannot return skewCorrectionVectors; mesh is now skewed"
            << abort(FatalError);
    }

    return (*skewCorrectionVectorsPtr_);
}


bool Foam::edgeInterpolation::movePoints() const
{
    lPNptr_.reset(nullptr);
    weightingFactorsPtr_.reset(nullptr);
    differenceFactorsPtr_.reset(nullptr);

    orthogonal_ = false;
    correctionVectorsPtr_.reset(nullptr);

    skew_ = true;
    skewCorrectionVectorsPtr_.reset(nullptr);

    return true;
}


const Foam::vector& Foam::edgeInterpolation::skewCorr(const label edgeI) const
{
    #ifdef FA_SKEW_CORRECTION

    return
        (
            skewCorrectionVectorsPtr_
          ? (*skewCorrectionVectorsPtr_)[edgeI]
          : pTraits<vector>::zero
        );

    #else

    return (*skewCorrectionVectorsPtr_)[edgeI];

    #endif
}


void Foam::edgeInterpolation::makeLPN() const
{
    DebugInFunction
        << "Constructing geodesic distance between points P and N"
        << endl;


    lPNptr_ = std::make_unique<surfaceScalarField>
    (
        IOobject
        (
            "lPN",
            mesh().time().constant(),
            mesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh(),
        dimLength
    );
    surfaceScalarField& lPN = *lPNptr_;

    // Set local references to mesh data
    const surfaceVectorField& edgeCentres = mesh().edgeCentres();
    const volVectorField& faceCentres = mesh().areaCentres();
    const labelUList& owner = mesh().owner();
    const labelUList& neighbour = mesh().neighbour();

    scalarField& lPNIn = lPN.primitiveFieldRef();

    // Calculate skewness correction vectors if necessary
    (void) skew();

    forAll(owner, edgeI)
    {
        const vector& skewCorrEdge = skewCorr(edgeI);

        scalar lPE =
            mag
            (
                edgeCentres[edgeI]
              - skewCorrEdge
              - faceCentres[owner[edgeI]]
            );

        scalar lEN =
            mag
            (
                faceCentres[neighbour[edgeI]]
              - edgeCentres[edgeI]
              + skewCorrEdge
            );

        lPNIn[edgeI] = (lPE + lEN);

        // Do not allow any mag(val) < SMALL
        if (mag(lPNIn[edgeI]) < SMALL)
        {
            lPNIn[edgeI] = SMALL;
        }
    }


    forAll(lPN.boundaryField(), patchI)
    {
        mesh().boundary()[patchI].makeDeltaCoeffs
        (
            lPN.boundaryFieldRef()[patchI]
        );

        lPN.boundaryFieldRef()[patchI] = 1.0/lPN.boundaryField()[patchI];
    }


    DebugInFunction
        << "Finished constructing geodesic distance PN"
        << endl;
}


void Foam::edgeInterpolation::makeWeights() const
{
    DebugInFunction
        << "Constructing weighting factors for edge interpolation"
        << endl;


    weightingFactorsPtr_ = std::make_unique<surfaceScalarField>
    (
        IOobject
        (
            "weightingFactors",
            mesh().pointsInstance(),
            mesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh(),
        dimensionedScalar(dimless, 1)
    );
    surfaceScalarField& weightingFactors = *weightingFactorsPtr_;


    // Set local references to mesh data
    const surfaceVectorField& edgeCentres = mesh().edgeCentres();
    const volVectorField& faceCentres = mesh().areaCentres();
    const labelUList& owner = mesh().owner();
    const labelUList& neighbour = mesh().neighbour();

    scalarField& weightingFactorsIn = weightingFactors.primitiveFieldRef();

    // Calculate skewness correction vectors if necessary
    (void) skew();

    forAll(owner, edgeI)
    {
        const vector& skewCorrEdge = skewCorr(edgeI);

        scalar lPE =
            mag
            (
                edgeCentres[edgeI]
              - skewCorrEdge
              - faceCentres[owner[edgeI]]
            );

        scalar lEN =
            mag
            (
                faceCentres[neighbour[edgeI]]
              - edgeCentres[edgeI]
              + skewCorrEdge
            );

        // weight = (0,1]
        const scalar lPN = lPE + lEN;
        if (mag(lPN) > SMALL)
        {
            weightingFactorsIn[edgeI] = lEN/lPN;
        }
    }

    forAll(mesh().boundary(), patchI)
    {
        mesh().boundary()[patchI].makeWeights
        (
            weightingFactors.boundaryFieldRef()[patchI]
        );
    }

    DebugInFunction
        << "Finished constructing weighting factors for face interpolation"
        << endl;
}


void Foam::edgeInterpolation::makeDeltaCoeffs() const
{
    DebugInFunction
        << "Constructing differencing factors array for edge gradient"
        << endl;

    // Force the construction of the weighting factors
    // needed to make sure deltaCoeffs are calculated for parallel runs.
    weights();

    differenceFactorsPtr_ = std::make_unique<surfaceScalarField>
    (
        IOobject
        (
            "differenceFactors",
            mesh().pointsInstance(),
            mesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh(),
        dimensionedScalar(dimless/dimLength, SMALL)
    );
    surfaceScalarField& DeltaCoeffs = *differenceFactorsPtr_;
    scalarField& dc = DeltaCoeffs.primitiveFieldRef();

    // Set local references to mesh data
    const surfaceVectorField& edgeCentres = mesh().edgeCentres();
    const volVectorField& faceCentres = mesh().areaCentres();
    const labelUList& owner = mesh().owner();
    const labelUList& neighbour = mesh().neighbour();
    const surfaceVectorField& lengths = mesh().Le();

//    const edgeList& edges = mesh().edges();
    const faceList& faces = mesh().faces();
    const pointField& points = mesh().points();

    // Calculate skewness correction vectors if necessary
    (void) skew();

    forAll(owner, edgeI)
    {
        const edge fe(faces[edgeI][0], faces[edgeI][1]);

        // Edge normal - area normal
        vector edgeNormal =
            normalised(lengths[edgeI] ^ fe.vec(points));

        Pout<< "makeDeltaCoeffs2 : edge:" << edgeI
            << " points:" << fe.line(points)
            << " vec:" << fe.vec(points)
            << " le:" << lengths[edgeI]
            << " normal:" << (lengths[edgeI] ^ fe.vec(points))
            << endl;

        // Unit delta vector
        vector unitDelta =
            faceCentres[neighbour[edgeI]]
          - faceCentres[owner[edgeI]];

        unitDelta.removeCollinear(edgeNormal);
        unitDelta.normalise();


        const vector& skewCorrEdge = skewCorr(edgeI);

        scalar lPE =
            mag
            (
                edgeCentres[edgeI]
              - skewCorrEdge
              - faceCentres[owner[edgeI]]
            );

        scalar lEN =
            mag
            (
                faceCentres[neighbour[edgeI]]
              - edgeCentres[edgeI]
              + skewCorrEdge
            );

        scalar lPN = lPE + lEN;


        // Edge normal - area tangent
        edgeNormal = normalised(lengths[edgeI]);

        // Do not allow any mag(val) < SMALL
        const scalar alpha = lPN*(unitDelta & edgeNormal);
        if (mag(alpha) > SMALL)
        {
            dc[edgeI] = scalar(1)/max(alpha, 0.05*lPN);
        }
    }


    forAll(DeltaCoeffs.boundaryField(), patchI)
    {
        mesh().boundary()[patchI].makeDeltaCoeffs
        (
            DeltaCoeffs.boundaryFieldRef()[patchI]
        );
    }
}


void Foam::edgeInterpolation::makeCorrectionVectors() const
{
    DebugInFunction
        << "Constructing non-orthogonal correction vectors"
        << endl;

    correctionVectorsPtr_ = std::make_unique<surfaceVectorField>
    (
        IOobject
        (
            "correctionVectors",
            mesh().pointsInstance(),
            mesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh(),
        dimless
    );
    surfaceVectorField& CorrVecs = *correctionVectorsPtr_;

    // Set local references to mesh data
    const volVectorField& faceCentres = mesh().areaCentres();

    const labelUList& owner = mesh().owner();
    const labelUList& neighbour = mesh().neighbour();

    const surfaceVectorField& lengths = mesh().Le();

//    const edgeList& edges = mesh().edges();
    const faceList& faces = mesh().faces();
    const pointField& points = mesh().points();

    scalarField deltaCoeffs(owner.size(), SMALL);

    vectorField& CorrVecsIn = CorrVecs.primitiveFieldRef();

    forAll(owner, edgeI)
    {
        const edge fe(faces[edgeI][0], faces[edgeI][1]);

        // Edge normal - area normal
        vector edgeNormal =
            normalised(lengths[edgeI] ^ fe.vec(points));

        Pout<< "makeCorrectionVectors2 : edge:" << edgeI
            << " points:" << fe.line(points)
            << " vec:" << fe.vec(points)
            << " le:" << lengths[edgeI]
            << " normal:" << (lengths[edgeI] ^ fe.vec(points))
            << endl;


        // Unit delta vector
        vector unitDelta =
            faceCentres[neighbour[edgeI]]
          - faceCentres[owner[edgeI]];

        unitDelta.removeCollinear(edgeNormal);
        unitDelta.normalise();

        // Edge normal - area tangent
        edgeNormal = normalised(lengths[edgeI]);

        // Do not allow any mag(val) < SMALL
        const scalar alpha = unitDelta & edgeNormal;
        if (mag(alpha) > SMALL)
        {
            deltaCoeffs[edgeI] = scalar(1)/alpha;
        }

        // Edge correction vector
        CorrVecsIn[edgeI] =
            edgeNormal
          - deltaCoeffs[edgeI]*unitDelta;
    }


    surfaceVectorField::Boundary& CorrVecsbf = CorrVecs.boundaryFieldRef();

    forAll(CorrVecs.boundaryField(), patchI)
    {
        //mesh().boundary()[patchI].makeCorrectionVectors(CorrVecsbf[patchI]);
        mesh().boundary()[patchI].makeNonOrthoCorrVectors(CorrVecsbf[patchI]);
    }

    scalar NonOrthogCoeff = 0.0;

    if (owner.size() > 0)
    {
        scalarField sinAlpha(deltaCoeffs*mag(CorrVecs.internalField()));
        sinAlpha.clamp_range(-1, 1);

        NonOrthogCoeff = max(Foam::asin(sinAlpha)*180.0/M_PI);
    }

    reduce(NonOrthogCoeff, maxOp<scalar>());

    DebugInFunction
        << "non-orthogonality coefficient = " << NonOrthogCoeff << " deg."
        << endl;

    if (NonOrthogCoeff < 0.1)
    {
        orthogonal_ = true;
        correctionVectorsPtr_.reset(nullptr);
    }
    else
    {
        orthogonal_ = false;
    }

    DebugInFunction
        << "Finished constructing non-orthogonal correction vectors"
        << endl;
}


void Foam::edgeInterpolation::makeSkewCorrectionVectors() const
{
    DebugInFunction
        << "Constructing skew correction vectors"
        << endl;

    skewCorrectionVectorsPtr_ = std::make_unique<surfaceVectorField>
    (
        IOobject
        (
            "skewCorrectionVectors",
            mesh().pointsInstance(),
            mesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh(),
        dimensionedVector(dimless, Zero)
    );
    surfaceVectorField& SkewCorrVecs = *skewCorrectionVectorsPtr_;

    // Set local references to mesh data
    const volVectorField& C = mesh().areaCentres();
    const surfaceVectorField& Ce = mesh().edgeCentres();

    const labelUList& owner = mesh().owner();
    const labelUList& neighbour = mesh().neighbour();

    const pointField& points = mesh().points();
    const faceList& faces = mesh().faces();


    forAll(neighbour, edgeI)
    {
        const vector& P = C[owner[edgeI]];
        const vector& N = C[neighbour[edgeI]];
        const edge fe(faces[edgeI][0], faces[edgeI][1]);
        const vector& S = points[fe.start()];

        // (T:Eq. 5.4)
        const vector d(N - P);
        const vector e(fe.vec(points));
        const vector de(d^e);
        const scalar alpha = magSqr(de);

        if (alpha < SMALL)
        {
            // Too small - skew correction remains zero
            continue;
        }
        const scalar beta = -((d^(S - P)) & de)/alpha;

        // (T:Eq. 5.3)
        const vector E(S + beta*e);

        SkewCorrVecs[edgeI] = Ce[edgeI] - E;
    }


    surfaceVectorField::Boundary& bSkewCorrVecs =
        SkewCorrVecs.boundaryFieldRef();

    forAll(SkewCorrVecs.boundaryField(), patchI)
    {
        fvsPatchVectorField& patchSkewCorrVecs = bSkewCorrVecs[patchI];

        if (patchSkewCorrVecs.coupled())
        {
            const labelUList& edgeFaces =
                mesh().boundary()[patchI].faceCells();

            //const edgeList::subList patchEdges =
            //    mesh().boundary()[patchI].patchSlice(edges);

            vectorField ngbC(C.boundaryField()[patchI].patchNeighbourField());

            forAll(patchSkewCorrVecs, edgeI)
            {
                const vector& P = C[edgeFaces[edgeI]];
                const vector& N = ngbC[edgeI];

                const face& f = mesh().boundary()[patchI].patch()[edgeI];
                const edge fe(f[0], f[1]);
                const vector& S = points[fe.start()];

                // (T:Eq. 5.4)
                const vector d(N - P);
                const vector e(fe.vec(points));
                const vector de(d^e);
                const scalar alpha = magSqr(de);

                if (alpha < SMALL)
                {
                    // Too small - skew correction remains zero
                    continue;
                }
                const scalar beta = -((d^(S - P)) & de)/alpha;

                const vector E(S + beta*e);

                patchSkewCorrVecs[edgeI] =
                    Ce.boundaryField()[patchI][edgeI] - E;
            }
        }
    }

    #ifdef FA_SKEW_CORRECTION

    constexpr scalar maxSkewRatio = 0.1;
    scalar skewCoeff = 0;

    forAll(own, edgeI)
    {
        const scalar magSkew = mag(skewCorrVecs[edgeI]);

        const scalar lPN =
            mag
            (
                Ce[edgeI]
              - skewCorrVecs[edgeI]
              - C[owner[edgeI]]
            )
          + mag
            (
                C[neighbour[edgeI]]
              - Ce[edgeI]
              + skewCorrVecs[edgeI]
            );

        const scalar ratio = magSkew/lPN;

        if (skewCoeff < ratio)
        {
            skewCoeff = ratio;

            if (skewCoeff > maxSkewRatio)
            {
                break;
            }
        }
    }

    DebugInFunction
        << "skew coefficient = " << skewCoeff << endl;

    if (skewCoeff < maxSkewRatio)
    {
        skewCorrectionVectorsPtr_.reset(nullptr);
    }

    #endif

    skew_ = bool(skewCorrectionVectorsPtr_);


    DebugInFunction
        << "Finished constructing skew correction vectors"
        << endl;
}


// ************************************************************************* //
