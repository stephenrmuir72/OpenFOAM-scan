/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Patrick Wilhelm, Data Lab Hell
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

#include "conformalFaceToFace.H"
#include "Time.H"
#include "faceSet.H"
#include "indexedOctree.H"
#include "treeDataPoint.H"
#include "mappedPatchBase.H"
#include "globalIndex.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(conformalFaceToFace, 0);
    addToRunTimeSelectionTable(topoSetSource, conformalFaceToFace, word);
    addToRunTimeSelectionTable(topoSetFaceSource, conformalFaceToFace, word);
}


Foam::topoSetSource::addToUsageTable Foam::conformalFaceToFace::usage_
(
    conformalFaceToFace::typeName,
    "\n    Usage: conformalFaceToFace region set\n\n"
    "    Select the boundary faces conformal with a faceSet in another region"
    "\n"
);


// Default tolerance as a fraction of source face size (centres) and area
static const Foam::scalar relTol = 1e-3;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::conformalFaceToFace::combine(topoSet& set, const bool add) const
{
    const polyMesh& meshSource = *meshSourcePtr_;
    const faceSet setSource(meshSource, sourceSetName_, IOobject::NO_REGISTER);

    const vectorField& CfSource = meshSource.faceCentres();
    const vectorField& SfSource = meshSource.faceAreas();

    // Source boundary faces on this rank: centre, unit normal and area
    DynamicList<point> localCf(setSource.size());
    DynamicList<vector> localNf(setSource.size());
    DynamicList<scalar> localArea(setSource.size());

    for (const label facei : setSource)
    {
        if
        (
            facei < meshSource.nInternalFaces()
         || facei >= meshSource.nFaces()
        )
        {
            FatalErrorInFunction
                << "Face " << facei << " in faceSet " << sourceSetName_
                << " is not a boundary face of region " << sourceRegionName_
                << exit(FatalError);
        }

        const scalar a = mag(SfSource[facei]);

        if (a < ROOTVSMALL)
        {
            FatalErrorInFunction
                << "Face " << facei << " in faceSet " << sourceSetName_
                << " has no area"
                << exit(FatalError);
        }

        localCf.append(CfSource[facei]);
        localNf.append(SfSource[facei]/a);
        localArea.append(a);
    }

    // Gather the source faces onto every rank, so each rank matches the whole
    // source against its own part of the target boundary
    List<point> CfSrc;
    List<vector> nfSrc;
    List<scalar> aSrc;
    globalIndex::gatherOp(localCf, CfSrc);
    globalIndex::gatherOp(localNf, nfSrc);
    globalIndex::gatherOp(localArea, aSrc);
    Pstream::broadcast(CfSrc);
    Pstream::broadcast(nfSrc);
    Pstream::broadcast(aSrc);

    const label nSrc = CfSrc.size();

    // Range of local target boundary faces to search
    const polyBoundaryMesh& patchesTarget = mesh_.boundaryMesh();

    label patchiTarget = -1;

    if (!patchName_.empty())
    {
        patchiTarget = patchesTarget.findPatchID(patchName_);

        if (patchiTarget < 0)
        {
            FatalErrorInFunction
                << "Cannot find patch " << patchName_
                << " in region " << mesh_.name() << nl
                << "Valid patches: " << patchesTarget.names()
                << exit(FatalError);
        }
    }

    const label faceiBegin =
    (
        patchiTarget < 0
      ? mesh_.nInternalFaces()
      : patchesTarget[patchiTarget].start()
    );
    const label faceiEnd =
    (
        patchiTarget < 0
      ? mesh_.nFaces()
      : faceiBegin + patchesTarget[patchiTarget].size()
    );

    const vectorField& CfTarget = mesh_.faceCentres();
    const vectorField& SfTarget = mesh_.faceAreas();

    pointField ctrsTarget(faceiEnd - faceiBegin);
    forAll(ctrsTarget, i)
    {
        ctrsTarget[i] = CfTarget[faceiBegin + i];
    }

    // A rank may hold no part of the target boundary; only the whole boundary
    // being empty is an error
    if (returnReduce(ctrsTarget.size(), sumOp<label>()) == 0)
    {
        FatalErrorInFunction
            << "No boundary faces to match against in region " << mesh_.name()
            << exit(FatalError);
    }

    const scalar distanceSigned = flip_ ? -distance_ : distance_;

    // The box must hold the largest per-face bound over all source faces
    scalar maxTol = tolerance_;
    if (!absTolerance_)
    {
        for (const scalar a : aSrc)
        {
            maxTol = max(maxTol, relTol*Foam::sqrt(a));
        }
    }

    const bool haveTarget = !ctrsTarget.empty();

    autoPtr<indexedOctree<treeDataPoint>> treePtr;
    if (haveTarget)
    {
        treeBoundBox bbTarget(ctrsTarget);
        bbTarget.grow(max(distance_, maxTol) + ROOTSMALL);

        treePtr.reset
        (
            new indexedOctree<treeDataPoint>
            (
                treeDataPoint(ctrsTarget),
                bbTarget,
                8,      // maxLevel
                10,     // leafsize
                3.0     // duplicity
            )
        );
    }

    // Per source face: the nearest matching target face on this rank and the
    // number of centres within the bound. A default nearInfo is a miss, which
    // loses to any hit under nearestEqOp
    List<mappedPatchBase::nearInfo> nearest(nSrc);
    labelList nNear(nSrc, Zero);

    if (haveTarget)
    {
        const indexedOctree<treeDataPoint>& treeTarget = treePtr();

        for (label g = 0; g < nSrc; ++g)
        {
            const point& Cf = CfSrc[g];
            const vector& nf = nfSrc[g];

            const scalar tol =
                absTolerance_
              ? tolerance_
              : relTol*Foam::sqrt(aSrc[g]);

            // Target search points: one signed separation along the normal
            // (directional), and one either way (anySign)
            const point pf = Cf + distanceSigned*nf;
            const point pfAlong = Cf + distance_*nf;
            const point pfAgainst = Cf - distance_*nf;

            scalar dMin = GREAT;
            label faceiMin = -1;

            if (directional_)
            {
                // Every candidate within the bound is weighed, so a second
                // face that close cannot mask the opposing one
                const labelList nearL = treeTarget.findSphere(pf, sqr(tol));
                nNear[g] = nearL.size();

                for (const label i : nearL)
                {
                    const label facej = faceiBegin + i;
                    const vector nt = normalised(SfTarget[facej]);
                    const scalar d = mag(CfTarget[facej] - pf);

                    // A conformal target opposes the source normal
                    if ((nt & nf) < 0 && d < dMin)
                    {
                        dMin = d;
                        faceiMin = facej;
                    }
                }
            }
            else
            {
                // Stacking is unknown, so the centre error is taken to the
                // nearer of the two points one separation along the normal
                const labelList nearL =
                    treeTarget.findSphere(Cf, sqr(distance_ + tol));
                nNear[g] = nearL.size();

                for (const label i : nearL)
                {
                    const label facej = faceiBegin + i;
                    const scalar d = min
                    (
                        mag(CfTarget[facej] - pfAlong),
                        mag(CfTarget[facej] - pfAgainst)
                    );

                    if (d < dMin)
                    {
                        dMin = d;
                        faceiMin = facej;
                    }
                }
            }

            if (faceiMin >= 0 && dMin <= tol)
            {
                nearest[g].first() =
                    pointIndexHit(true, CfTarget[faceiMin], faceiMin);
                nearest[g].second() =
                    Tuple2<scalar, label>(sqr(dMin), UPstream::myProcNo());
            }
        }
    }

    Pstream::listCombineReduce(nearest, mappedPatchBase::nearestEqOp());
    Pstream::listReduce(nNear, sumOp<label>());

    // Every source face must match a target face somewhere
    forAll(nearest, g)
    {
        if (nearest[g].first().hit())
        {
            continue;
        }

        const scalar tol =
            absTolerance_ ? tolerance_ : relTol*Foam::sqrt(aSrc[g]);

        FatalErrorInFunction
            << "No target face matches source face at " << CfSrc[g]
            << " within tolerance " << tol << nl;

        if (directional_ && nNear[g])
        {
            FatalError
                << nNear[g] << " target centres lie within the bound but none"
                << " opposes the source face normal " << nfSrc[g];
        }

        FatalError << exit(FatalError);
    }

    // Add the matched target faces on their owning rank, checking the area and
    // one-to-one matching there
    labelHashSet matched(2*nSrc);
    label nAreaMismatch = 0;
    scalar areaSource = 0;
    scalar areaTarget = 0;

    forAll(nearest, g)
    {
        if (nearest[g].second().second() != UPstream::myProcNo())
        {
            continue;
        }

        // A centre match alone accepts a refined sub-face at the same centre
        const label faceiMin = nearest[g].first().index();
        const scalar aSource = aSrc[g];
        const scalar aTarget = mag(SfTarget[faceiMin]);
        const scalar aDiff = mag(aTarget - aSource)/aSource;

        if (aDiff > areaTolerance_)
        {
            FatalErrorInFunction
                << "Source face with centre " << CfSrc[g]
                << " matched target face " << faceiMin << nl
                << "The areas differ by " << aDiff << " of the source area "
                << aSource << "; areaTolerance bounds the difference"
                << exit(FatalError);
        }

        if (aDiff > relTol)
        {
            ++nAreaMismatch;
        }

        areaSource += aSource;
        areaTarget += aTarget;

        if (!matched.insert(faceiMin))
        {
            FatalErrorInFunction
                << "Target face " << faceiMin
                << " matches more than one source face"
                << exit(FatalError);
        }

        addOrDelete(set, faceiMin, add);
    }

    reduce(nAreaMismatch, sumOp<label>());
    reduce(areaSource, sumOp<scalar>());
    reduce(areaTarget, sumOp<scalar>());

    if (nAreaMismatch)
    {
        WarningInFunction
            << nAreaMismatch << " of " << nSrc
            << " matched pairs differ in area" << nl
            << "The matched target faces carry " << areaTarget
            << " of the " << areaSource << " carried by the source faces"
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::conformalFaceToFace::conformalFaceToFace
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetFaceSource(mesh, dict),
    sourceRegionName_(dict.get<word>("region")),
    sourceSetName_(dict.get<word>("set")),
    patchName_(dict.getOrDefault<word>("patch", word::null)),
    distance_
    (
        dict.getCheckOrDefault<scalar>("distance", 0, scalarMinMax::ge(0))
    ),
    absTolerance_(dict.found("tolerance")),
    tolerance_
    (
        dict.getCheckOrDefault<scalar>("tolerance", 0, scalarMinMax::ge(0))
    ),
    areaTolerance_
    (
        dict.getCheckOrDefault<scalar>
        (
            "areaTolerance", relTol, scalarMinMax::ge(0)
        )
    ),
    directional_(!dict.getOrDefault("anySign", false)),
    flip_(dict.getOrDefault("flipDistance", false)),
    meshSourcePtr_(nullptr)
{
    if (!directional_ && flip_)
    {
        FatalErrorInFunction
            << "flipDistance sets the search direction, which anySign discards"
            << exit(FatalError);
    }

    if (flip_ && distance_ == 0)
    {
        FatalErrorInFunction
            << "flipDistance has no separation to reverse at zero distance"
            << exit(FatalError);
    }

    meshSourcePtr_.reset
    (
        new polyMesh
        (
            IOobject
            (
                sourceRegionName_,
                mesh.time().timeName(),
                mesh.time(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        )
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::conformalFaceToFace::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding boundary faces conformal with faceSet "
                << sourceSetName_ << " of region " << sourceRegionName_
                << " ..." << endl;
        }

        combine(set, true);
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing boundary faces conformal with faceSet "
                << sourceSetName_ << " of region " << sourceRegionName_
                << " ..." << endl;
        }

        combine(set, false);
    }
}


// ************************************************************************* //
