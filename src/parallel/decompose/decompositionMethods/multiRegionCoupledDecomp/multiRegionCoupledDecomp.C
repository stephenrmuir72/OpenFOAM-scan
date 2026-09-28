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

#include "multiRegionCoupledDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "regionProperties.H"
#include "mappedPatchBase.H"
#include "polyMesh.H"
#include "Time.H"
#include "IOdictionary.H"
#include "DynamicList.H"
#include "labelPair.H"
#include "ListOps.H"
#include "MinMax.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiRegionCoupledDecomp, 0);
    addToRunTimeSelectionTable
    (
        decompositionMethod,
        multiRegionCoupledDecomp,
        dictionary
    );

    defineTypeName(multiRegionCoupledDecomp::weldedDecomposition);
}


Foam::multiRegionCoupledDecomp::weldedDecomposition::weldedDecomposition
(
    const objectRegistry& obr
)
:
    regIOobject
    (
        IOobject
        (
            typeName,
            obr.time().timeName(),
            obr,
            IOobjectOption::NO_READ,
            IOobjectOption::NO_WRITE,
            IOobjectOption::REGISTER
        )
    )
{}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Union-find with path halving; smaller root index wins for determinism
static inline label findRoot(labelUList& parent, label celli)
{
    while (parent[celli] != celli)
    {
        parent[celli] = parent[parent[celli]];
        celli = parent[celli];
    }
    return celli;
}

static inline void unionCells(labelUList& parent, const label a, const label b)
{
    const label roota = findRoot(parent, a);
    const label rootb = findRoot(parent, b);
    if (roota < rootb)
    {
        parent[rootb] = roota;
    }
    else if (rootb < roota)
    {
        parent[roota] = rootb;
    }
}

// The dictionary decompositionModel reads for a region,
// system/<region>/decomposeParDict
static IOobject regionDecompDictIO(const polyMesh& mesh)
{
    return IOobject
    (
        "decomposeParDict",     // decompositionModel::canonicalName
        mesh.time().system(),
        mesh.local(),
        mesh.thisDb(),
        IOobjectOption::MUST_READ,
        IOobjectOption::NO_WRITE,
        IOobjectOption::NO_REGISTER,
        true    // global object
    );
}


// Compare dictionaries by content, independent of the entry order
static bool sameDict(const dictionary& a, const dictionary& b)
{
    const wordList keys(a.sortedToc());

    if (keys != b.sortedToc())
    {
        return false;
    }

    for (const word& key : keys)
    {
        const entry& ea = *a.findEntry(key, keyType::LITERAL);
        const entry& eb = *b.findEntry(key, keyType::LITERAL);

        if (ea.isDict() != eb.isDict())
        {
            return false;
        }

        if (ea.isDict())
        {
            if (!sameDict(ea.dict(), eb.dict()))
            {
                return false;
            }
        }
        else
        {
            const tokenList& ta = ea.stream();
            const tokenList& tb = eb.stream();

            if (ta != tb)
            {
                return false;
            }
        }
    }

    return true;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiRegionCoupledDecomp::multiRegionCoupledDecomp
(
    const dictionary& decompDict,
    const word& regionName
)
:
    decompositionMethod(decompDict, regionName),
    methodDict_
    (
        findCoeffsDict
        (
            typeName + "Coeffs",
            (selectionType::EXACT | selectionType::MANDATORY)
        )
    ),
    method_(nullptr)
{
    if (nDomains() < 1)
    {
        FatalIOErrorInFunction(decompDict_)
            << "numberOfSubdomains is " << nDomains()
            << "; multiRegionCoupled sizes its per-processor lists from"
            << " numberOfSubdomains, which must be at least 1"
            << exit(FatalIOError);
    }

    if (!constraints_.empty())
    {
        FatalIOErrorInFunction(decompDict_)
            << typeName << " does not support decomposition constraints:"
            << " applied per region, they can split a welded pair across"
            << " processors"
            << exit(FatalIOError);
    }

    methodDict_.set("numberOfSubdomains", nDomains());
    method_ = decompositionMethod::New(methodDict_);
    methodDict_.readIfPresent("implicitPatches", implicitPatches_);
    methodDict_.readIfPresent("explicitPatches", explicitPatches_);
    if
    (
        const auto* eptr =
            methodDict_.findEntry("regionWeights", keyType::LITERAL)
    )
    {
        if (const auto* dictptr = eptr->dictPtr())
        {
            regionWeights_ = *dictptr;
        }
        else
        {
            FatalIOErrorInFunction(methodDict_)
                << "regionWeights is not a dictionary" << exit(FatalIOError);
        }
    }

    // The regions entry of this file; a region's own decomposeParDict is
    // checked in checkRegionMethods, once the region list is known
    if (const auto* regionsPtr = decompDict_.findDict("regions"))
    {
        for (const entry& dEntry : *regionsPtr)
        {
            word regionMethod;
            const auto* dictptr = dEntry.dictPtr();
            if
            (
                dictptr
             && dictptr->readIfPresent("method", regionMethod)
             && regionMethod != typeName
            )
            {
                FatalIOErrorInFunction(*regionsPtr)
                    << "Region " << dEntry.keyword() << " overrides method"
                    << " with " << regionMethod << "; multiRegionCoupled"
                    << " balances the cells of every region in"
                    << " regionProperties as one decomposition, so no region"
                    << " can be placed by another method"
                    << exit(FatalIOError);
            }
        }
    }
    if
    (
        const auto* eptr =
            methodDict_.findEntry("patchWeights", keyType::LITERAL)
    )
    {
        const auto* dictptr = eptr->dictPtr();
        if (!dictptr)
        {
            FatalIOErrorInFunction(methodDict_)
                << "patchWeights is not a dictionary" << exit(FatalIOError);
        }

        patchWeights_ = *dictptr;

        for (const entry& dEntry : patchWeights_)
        {
            if (!dEntry.dictPtr())
            {
                FatalIOErrorInFunction(patchWeights_)
                    << "patchWeights entry " << dEntry.keyword()
                    << " is not a sub-dictionary; the boundary-face weights "
                    << "are given per region" << exit(FatalIOError);
            }
        }
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::multiRegionCoupledDecomp::weldedDecomposition&
Foam::multiRegionCoupledDecomp::weldedDecomp(const polyMesh& mesh) const
{
    const objectRegistry& obr = mesh.time();

    auto* weldPtr =
        obr.getObjectPtr<weldedDecomposition>(weldedDecomposition::typeName);

    if (!weldPtr)
    {
        weldPtr = new weldedDecomposition(obr);
        weldPtr->ownerName = mesh.name();
        weldPtr->coeffs = methodDict_;
        regIOobject::store(weldPtr);
        buildGlobalDecomposition(mesh, *weldPtr);
    }
    // numberOfSubdomains is injected into methodDict_, so sameDict below
    // would catch a mismatch too; check it first for a message that names
    // the counts rather than pointing at two identical-looking Coeffs
    else if
    (
        const label nWeldDomains =
            weldPtr->coeffs.get<label>("numberOfSubdomains");
        nWeldDomains != nDomains()
    )
    {
        FatalErrorInFunction
            << "Region " << mesh.name() << " asks for " << nDomains()
            << " subdomains, region " << weldPtr->ownerName << " built the"
            << " decomposition for " << nWeldDomains
            << "; one decomposition balances every region in regionProperties"
            << exit(FatalError);
    }
    else if (!sameDict(methodDict_, weldPtr->coeffs))
    {
        // decomposePar resolves the coefficients per region, from
        // system/<region>/decomposeParDict or a regions sub-dictionary,
        // but one decomposition balances every region in regionProperties
        FatalErrorInFunction
            << "Region " << mesh.name() << " gives different " << typeName
            << "Coeffs than region " << weldPtr->ownerName << ", which built "
            << "the decomposition; every region must give the same "
            << "coefficients" << exit(FatalError);
    }

    return *weldPtr;
}


void Foam::multiRegionCoupledDecomp::buildGlobalDecomposition
(
    const polyMesh& mesh,
    weldedDecomposition& weld
) const
{
    const Time& runTime = mesh.time();

    wordList& regions = weld.regions;
    labelList& regionOffsets = weld.regionOffsets;
    labelList& globalDecomp = weld.decomp;

    // All regions, in regionProperties order
    regions = regionProperties(runTime).names();
    const label nRegions = regions.size();

    HashTable<label> regionId;
    forAll(regions, regioni)
    {
        regionId.insert(regions[regioni], regioni);
    }

    // Current region is already loaded; load the others read-only and
    // unregistered, so the mapped sampleMesh lookup is never triggered.
    PtrList<polyMesh> ownMeshes(nRegions);
    UPtrList<const polyMesh> meshes(nRegions);
    forAll(regions, regioni)
    {
        if (regions[regioni] == mesh.name())
        {
            meshes.set(regioni, &mesh);
        }
        else
        {
            ownMeshes.set
            (
                regioni,
                new polyMesh
                (
                    IOobject
                    (
                        regions[regioni],
                        runTime.timeName(),
                        runTime,
                        IOobjectOption::MUST_READ,
                        IOobjectOption::NO_WRITE,
                        IOobjectOption::NO_REGISTER
                    )
                )
            );
            meshes.set(regioni, ownMeshes.get(regioni));
        }
    }

    checkRegionMethods(mesh, meshes);

    // Global cell offsets (cf. lduPrimitiveMeshAssembly cellOffsets_)
    regionOffsets.setSize(nRegions + 1);
    regionOffsets[0] = 0;
    forAll(regions, regioni)
    {
        regionOffsets[regioni+1] =
            regionOffsets[regioni] + meshes[regioni].nCells();
    }
    const label nGlobal = regionOffsets[nRegions];

    // Weld the cells either side of every coupled face; the pairs kept out
    // of the weld are recorded to check the decomposition below. Each pair
    // is seen from both sides; union is idempotent.
    labelList parent(identity(nGlobal));
    DynamicList<labelPair> excluded;
    label nWelded = 0;

    const wordRes::filter weldFilter(implicitPatches_, explicitPatches_);

    forAll(regions, regioni)
    {
        const polyBoundaryMesh& pbm = meshes[regioni].boundaryMesh();

        forAll(pbm, patchi)
        {
            const polyPatch& pp = pbm[patchi];

            const auto* mppPtr = isA<mappedPatchBase>(pp);
            if (!mppPtr)
            {
                continue;
            }

            const mappedPatchBase& mpp = *mppPtr;

            // A patch the filters exclude by its own name is skipped before
            // samplePatch() is read: a nearestCell patch carries none and
            // the accessor aborts on it
            const bool ownExcluded = !weldFilter(pp.name());

            if (!mpp.coupleGroup().empty())
            {
                if (ownExcluded)
                {
                    continue;
                }

                FatalErrorInFunction
                    << "Coupled patch " << regions[regioni] << '/'
                    << pp.name() << " is wired by coupleGroup; specify"
                    << " sampleRegion and samplePatch" << exit(FatalError);
            }

            const label nbrRegioni = regionId.lookup(mpp.sampleRegion(), -1);
            if (nbrRegioni < 0)
            {
                continue;   // samples a region outside the coupled set
            }

            // The assembly pairs face i with neighbour face i, so only a
            // face-to-face sample mode gives a meaningful pairing
            if (mpp.mode() != mappedPatchBase::NEARESTPATCHFACE)
            {
                if (ownExcluded)
                {
                    continue;
                }

                FatalErrorInFunction
                    << "Coupled patch " << regions[regioni] << '/' << pp.name()
                    << " is not nearestPatchFace; implicit multi-region "
                    << "coupling pairs face i with neighbour face i"
                    << exit(FatalError);
            }

            const word& nbrName = mpp.samplePatch();

            // weldFilter over both sides of the pair
            const auto ownMatch = wordRes::matched(implicitPatches_, pp.name());
            const auto nbrMatch = wordRes::matched(implicitPatches_, nbrName);

            const bool welded =
                ownMatch == wordRe::LITERAL
             || nbrMatch == wordRe::LITERAL
             || (
                    (
                        implicitPatches_.empty()
                     || ownMatch == wordRe::REGEX
                     || nbrMatch == wordRe::REGEX
                    )
                 && !explicitPatches_.match(pp.name())
                 && !explicitPatches_.match(nbrName)
                );

            const label nbrPatchi =
                meshes[nbrRegioni].boundaryMesh().findPatchID(nbrName);
            if (nbrPatchi < 0)
            {
                if (!welded)
                {
                    continue;
                }

                FatalErrorInFunction
                    << "Sample patch " << mpp.sampleRegion() << '/' << nbrName
                    << " for coupled patch " << regions[regioni] << '/'
                    << pp.name() << " not found" << exit(FatalError);
            }

            const labelUList& fc = pp.faceCells();
            const labelUList& nbrFc =
                meshes[nbrRegioni].boundaryMesh()[nbrPatchi].faceCells();

            if (fc.size() != nbrFc.size())
            {
                if (!welded)
                {
                    continue;
                }

                FatalErrorInFunction
                    << "Coupled patches " << regions[regioni] << '/'
                    << pp.name() << " (" << fc.size() << " faces) and "
                    << mpp.sampleRegion() << '/' << nbrName
                    << " (" << nbrFc.size() << " faces) are not conformal"
                    << exit(FatalError);
            }

            // Face i to neighbour face i, as mapCollocatedFaces() pairs them;
            // the centres are not checked. An interface whose two sides map
            // back at each other is handled from the lower (region, patch).
            const polyPatch& nbrPp =
                meshes[nbrRegioni].boundaryMesh()[nbrPatchi];

            bool twoSided = false;
            if (const auto* nbrMppPtr = isA<mappedPatchBase>(nbrPp))
            {
                twoSided =
                    nbrMppPtr->mode() == mappedPatchBase::NEARESTPATCHFACE
                 && nbrMppPtr->coupleGroup().empty()
                 && nbrMppPtr->sampleRegion() == regions[regioni]
                 && nbrMppPtr->samplePatch() == pp.name();
            }

            const bool onceOnly =
                !twoSided
             || regioni < nbrRegioni
             || (regioni == nbrRegioni && patchi <= nbrPatchi);

            if (welded && onceOnly)
            {
                nWelded += fc.size();
            }

            forAll(fc, i)
            {
                const label globali = regionOffsets[regioni] + fc[i];
                const label nbrGlobali = regionOffsets[nbrRegioni] + nbrFc[i];

                if (welded)
                {
                    unionCells(parent, globali, nbrGlobali);
                }
                else if (onceOnly)
                {
                    excluded.emplace_back(globali, nbrGlobali);
                }
            }
        }
    }

    // Contract welded cells to dense super-cells
    labelList superOf(nGlobal, -1);
    label nSuper = 0;
    forAll(superOf, globali)
    {
        const label root = findRoot(parent, globali);
        if (superOf[root] < 0)
        {
            superOf[root] = nSuper++;
        }
    }
    forAll(superOf, globali)
    {
        superOf[globali] = superOf[findRoot(parent, globali)];
    }

    // Per-region cost factors (default 1) scale the per-cell cost
    scalarList regionWeight(nRegions, scalar(1));
    forAll(regions, regioni)
    {
        regionWeights_.readIfPresent(regions[regioni], regionWeight[regioni]);

        if (regionWeight[regioni] <= 0)
        {
            FatalIOErrorInFunction(regionWeights_)
                << "regionWeights entry for " << regions[regioni] << " is "
                << regionWeight[regioni] << "; a weight has to be positive"
                << exit(FatalIOError);
        }
    }

    // readIfPresent matches patterns, so the membership test does too
    DynamicList<word> unmatched(regionWeights_.size());
    for (const entry& dEntry : regionWeights_)
    {
        const keyType& key = dEntry.keyword();

        bool matched = false;
        forAll(regions, regioni)
        {
            if (key.match(regions[regioni]))
            {
                matched = true;
                break;
            }
        }

        if (!matched)
        {
            unmatched.push_back(key);
        }
    }
    if (!unmatched.empty())
    {
        WarningInFunction
            << "regionWeights entries " << unmatched
            << " match no region" << endl;
    }

    // Per-face cost factors of the region's patches, charged to the cell
    // owning the face
    const bool patchWeighted = !patchWeights_.empty();
    scalarField cellFaceLoad;
    labelList cellFaceCount;
    labelList regionFaces(nRegions, Zero);
    label nWeightedFaces = 0;

    if (patchWeighted)
    {
        cellFaceLoad.resize_fill(nGlobal, Zero);
        cellFaceCount.resize_fill(nGlobal, Zero);

        forAll(regions, regioni)
        {
            const auto* patchDictPtr =
                patchWeights_.findDict(regions[regioni]);

            if (!patchDictPtr)
            {
                continue;
            }

            for (const polyPatch& pp : meshes[regioni].boundaryMesh())
            {
                scalar w = 0;
                if (!patchDictPtr->readIfPresent(pp.name(), w))
                {
                    continue;
                }

                if (w < 0)
                {
                    FatalIOErrorInFunction(*patchDictPtr)
                        << "Negative weight " << w << " for patch "
                        << regions[regioni] << '/' << pp.name()
                        << exit(FatalIOError);
                }

                regionFaces[regioni] += pp.size();
                nWeightedFaces += pp.size();

                for (const label celli : pp.faceCells())
                {
                    const label globali = regionOffsets[regioni] + celli;
                    cellFaceLoad[globali] += w;
                    ++cellFaceCount[globali];
                }
            }
        }
    }

    const bool weighted = !regionWeights_.empty() || patchWeighted;

    // Super-cell weights = summed cost of the welded member cells
    scalarField superWeights(nSuper, Zero);
    forAll(regions, regioni)
    {
        const scalar w = regionWeight[regioni];
        const label start = regionOffsets[regioni];
        const label end = regionOffsets[regioni+1];
        for (label globali = start; globali < end; ++globali)
        {
            superWeights[superOf[globali]] += w;
        }
    }
    forAll(cellFaceLoad, globali)
    {
        superWeights[superOf[globali]] += cellFaceLoad[globali];
    }

    // Super-cell adjacency from the intra-region cell connectivity, counted
    // then filled, since a hash set per super-cell costs a pointer table and
    // a heap node per edge
    labelListList superCellCells(nSuper);
    {
        labelList nAdj(nSuper, Zero);
        forAll(regions, regioni)
        {
            const labelListList& cellCells = meshes[regioni].cellCells();
            forAll(cellCells, celli)
            {
                nAdj[superOf[regionOffsets[regioni] + celli]] +=
                    cellCells[celli].size();
            }
        }

        forAll(superCellCells, superi)
        {
            superCellCells[superi].resize_nocopy(nAdj[superi]);
        }

        nAdj = Zero;
        forAll(regions, regioni)
        {
            const labelListList& cellCells = meshes[regioni].cellCells();
            forAll(cellCells, celli)
            {
                const label superi = superOf[regionOffsets[regioni] + celli];
                for (const label nbrCelli : cellCells[celli])
                {
                    const label nbrSuperi =
                        superOf[regionOffsets[regioni] + nbrCelli];
                    if (superi != nbrSuperi)
                    {
                        superCellCells[superi][nAdj[superi]++] = nbrSuperi;
                    }
                }
            }
        }

        // Drop neighbours reached through more than one member cell; the
        // sort also keeps the row independent of the traversal order
        forAll(superCellCells, superi)
        {
            superCellCells[superi].resize(nAdj[superi]);
            inplaceUniqueSort(superCellCells[superi]);
        }
    }

    // Centroids for geometric sub-methods (graph methods ignore them)
    pointField superCtr(nSuper, Zero);
    labelList superCount(nSuper, Zero);
    forAll(regions, regioni)
    {
        const pointField& cc = meshes[regioni].cellCentres();
        forAll(cc, celli)
        {
            const label superi = superOf[regionOffsets[regioni] + celli];
            superCtr[superi] += cc[celli];
            ++superCount[superi];
        }
    }
    forAll(superCtr, superi)
    {
        if (superCount[superi])
        {
            superCtr[superi] /= superCount[superi];
        }
    }

    // Load-balance the welded super-cells, then expand to per-cell processors
    const labelList superProc
    (
        method_().decompose(superCellCells, superCtr, superWeights)
    );

    globalDecomp.setSize(nGlobal);
    forAll(globalDecomp, globali)
    {
        globalDecomp[globali] = superProc[superOf[globali]];
    }

    // Both cells of a coupled face must share a processor for an implicit
    // (coupledEnergyField) run. lduPrimitiveMeshAssembly sees the split as
    // max(iFaces - nbrFaces, 0) on the master side, so it catches only the
    // ranks where the master patch is the larger one.
    label nSplit = 0;
    for (const labelPair& facePair : excluded)
    {
        if (globalDecomp[facePair.first()] != globalDecomp[facePair.second()])
        {
            ++nSplit;
        }
    }
    if (nSplit > 0)
    {
        WarningInFunction
            << nSplit << '/' << excluded.size() << " coupled faces left out"
            << " of the weld have their cells on different processors; an"
            << " implicit (coupledEnergyField) run over those patches either"
            << " aborts in lduPrimitiveMeshAssembly or pairs the two sides by"
            << " patch-local index" << endl;
    }

    if (!implicitPatches_.empty())
    {
        Info<< type() << " : welding only patches matching "
            << implicitPatches_ << endl;
    }
    if (!explicitPatches_.empty())
    {
        Info<< type() << " : excluding patches matching "
            << explicitPatches_ << endl;
    }
    if (patchWeighted)
    {
        if (nWeightedFaces)
        {
            Info<< type() << " : weighting " << nWeightedFaces
                << " boundary faces (";
            forAll(regions, regioni)
            {
                if (regioni)
                {
                    Info<< ' ';
                }
                Info<< regions[regioni] << ':' << regionFaces[regioni];
            }
            Info<< ')' << endl;
        }
        else
        {
            WarningInFunction
                << "patchWeights matched no boundary face" << endl;
        }
    }

    // Largest indivisible super-cell; bounds the achievable balance
    const label maxSuper = max(superCount);

    Info<< type() << " : welded " << nGlobal << " cells across "
        << nRegions << " regions into " << nSuper << " super-cells (largest "
        << maxSuper << " cells, " << nWelded
        << " coupled faces welded)" << endl;

    // Global cross-region load balance; not shown by the per-region summaries
    const label nProcs = nDomains();

    labelList procCells(nProcs, Zero);
    labelList procFaces(nProcs, Zero);
    scalarField procLoad(nProcs, Zero);
    labelListList procRegionCells(nProcs, labelList(nRegions, Zero));
    forAll(regions, regioni)
    {
        const label start = regionOffsets[regioni];
        const label end = regionOffsets[regioni+1];
        for (label globali = start; globali < end; ++globali)
        {
            const label proci = globalDecomp[globali];
            ++procCells[proci];
            ++procRegionCells[proci][regioni];
            procLoad[proci] += regionWeight[regioni];
        }
    }
    forAll(cellFaceCount, globali)
    {
        const label proci = globalDecomp[globali];
        procFaces[proci] += cellFaceCount[globali];
        procLoad[proci] += cellFaceLoad[globali];
    }

    {
        DynamicList<labelPair> starved;
        forAll(procRegionCells, proci)
        {
            forAll(regions, regioni)
            {
                if (!procRegionCells[proci][regioni])
                {
                    starved.emplace_back(proci, regioni);
                }
            }
        }

        if (!starved.empty())
        {
            OSstream& os = WarningInFunction;
            os  << starved.size() << " (processor, region) pair(s) hold no"
                << " cells:" << nl;

            for (const labelPair& pr : starved)
            {
                os  << "    processor " << pr.first() << ": "
                    << regions[pr.second()] << nl;
            }

            os  << endl;
        }
    }

    const MinMax<label> cellRange = minMax(procCells);
    const scalar meanCells = scalar(nGlobal)/nProcs;

    Info<< type() << " : global load balance over " << nProcs
        << " processors" << nl;
    forAll(procCells, proci)
    {
        Info<< "    Processor " << proci << "  cells:" << procCells[proci]
            << " (";
        forAll(regions, regioni)
        {
            if (regioni)
            {
                Info<< ' ';
            }
            Info<< regions[regioni] << ':' << procRegionCells[proci][regioni];
        }
        Info<< ')';
        if (patchWeighted)
        {
            Info<< "  faces:" << procFaces[proci];
        }
        Info<< nl;
    }
    Info<< "    Cells min:" << cellRange.min() << " max:" << cellRange.max()
        << " mean:" << meanCells;
    if (meanCells > 0)
    {
        Info<< " imbalance:" << cellRange.max()/meanCells - 1;
    }

    if (patchWeighted)
    {
        const MinMax<label> faceRange = minMax(procFaces);
        const scalar meanFaces = scalar(nWeightedFaces)/nProcs;

        Info<< nl << "    Faces min:" << faceRange.min()
            << " max:" << faceRange.max() << " mean:" << meanFaces;
        if (meanFaces > 0)
        {
            Info<< " imbalance:" << faceRange.max()/meanFaces - 1;
        }
    }

    if (weighted)
    {
        const MinMax<scalar> loadRange = minMax(procLoad);
        const scalar meanLoad = sum(procLoad)/nProcs;

        Info<< nl << "    Load  min:" << loadRange.min()
            << " max:" << loadRange.max() << " mean:" << meanLoad;
        if (meanLoad > 0)
        {
            Info<< " imbalance:" << loadRange.max()/meanLoad - 1;
        }
    }
    Info<< endl;
}


void Foam::multiRegionCoupledDecomp::checkRegionMethods
(
    const polyMesh& mesh,
    const UPtrList<const polyMesh>& meshes
) const
{
    // A -decomposeParDict override hands every region the same file, which
    // the regions guard of the constructor covers
    if (decompDict_.name() != regionDecompDictIO(mesh).objectPath())
    {
        return;
    }

    for (const polyMesh& regionMesh : meshes)
    {
        if (regionMesh.name() == mesh.name())
        {
            continue;
        }

        const dictionary regionDict
        (
            IOdictionary::readContents(regionDecompDictIO(regionMesh))
        );

        // decompositionMethod::New lets a regions sub-dictionary of that
        // same file override the top-level method
        word regionMethod(regionDict.get<word>("method"));

        if
        (
            const auto* regionsPtr =
                cfindRegionDict(regionDict, regionMesh.name())
        )
        {
            regionsPtr->readIfPresent("method", regionMethod);
        }

        if (regionMethod != typeName)
        {
            FatalIOErrorInFunction(regionDict)
                << "Region " << regionMesh.name() << " is decomposed by "
                << regionMethod << "; multiRegionCoupled balances the cells"
                << " of every region in regionProperties as one"
                << " decomposition, so no region can be placed by another"
                << " method"
                << exit(FatalIOError);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::multiRegionCoupledDecomp::decompose
(
    const polyMesh& mesh,
    const pointField& points,
    const scalarField& cWeights
) const
{
    if (UPstream::parRun())
    {
        FatalErrorInFunction
            << typeName << " reads every region mesh and numbers the cells"
            << " globally across regions, so it runs in serial only"
            << exit(FatalError);
    }

    if (!cWeights.empty())
    {
        FatalErrorInFunction
            << "weightField is not used by " << typeName
            << ": one decomposition balances every region in regionProperties,"
            << " weighted through the coefficients of this method"
            << exit(FatalError);
    }

    const weldedDecomposition& weld = weldedDecomp(mesh);

    const label regioni = weld.regions.find(mesh.name());
    if (regioni < 0)
    {
        FatalErrorInFunction
            << "Region " << mesh.name()
            << " is not listed in constant/regionProperties"
            << exit(FatalError);
    }

    const label start = weld.regionOffsets[regioni];
    const label nCells = weld.regionOffsets[regioni+1] - start;

    if (nCells != points.size())
    {
        FatalErrorInFunction
            << "Cell count for region " << mesh.name()
            << " changed since the coupled decomposition was built: "
            << nCells << " vs " << points.size()
            << exit(FatalError);
    }

    labelList finalDecomp(nCells);
    forAll(finalDecomp, celli)
    {
        finalDecomp[celli] = weld.decomp[start + celli];
    }

    return finalDecomp;
}


// ************************************************************************* //
