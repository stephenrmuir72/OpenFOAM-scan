/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
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

#include "RugeStubenCoarsening.H"
#include "Map.H"
#include "DynamicList.H"
#include "UPstream.H"
#include "UOPstream.H"
#include "UIPstream.H"
#include "PstreamReduceOps.H"
#include "processorLduInterface.H"
#include <algorithm>
#include <queue>
#include <utility>
#include <vector>

#include <map>

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * //

namespace
{

// Exchange a per-interface labelList (already sized to each interface's
// own face count, in that interface's own face order) directly between
// the two ranks of every processor interface, using raw two-sided
// nonBlocking read/write. The data is already per-face (not per-cell),
// so no faceCells extraction is needed - unlike
// lduInterface::initInternalFieldTransfer/internalFieldTransfer, which
// are per-cell-field oriented.
//
// NOTE on comms pattern: processorLduInterface::send() (nonBlocking) bundles
// an internal receive-post sized to the LOCAL send buffer's OWN byte size
// (it assumes the incoming message is the same size as the outgoing one).
// That is only safe for perfectly symmetric, fixed-size, same-size-on-both-
// -sides exchanges. To avoid any risk of an MPI_ERR_TRUNCATE (and to match
// the exact idiom used elsewhere, e.g.
// processorGAMGInterfaceField::initInterfaceMatrixUpdate()), all exchanges
// here instead go through the raw UIPstream::read()/UOPstream::write()
// statics directly, posting every receive (sized explicitly and
// correctly by the caller) BEFORE posting any send.
Foam::List<Foam::labelList> exchangePerFaceLabel
(
    const Foam::lduInterfaceFieldPtrsList& interfaces,
    const Foam::List<Foam::labelList>& localPerFace
)
{
    using namespace Foam;

    const label startOfRequests = UPstream::nRequests();

    List<labelList> nbrPerFace(interfaces.size());

    forAll(interfaces, patchi)
    {
        if
        (
            interfaces.set(patchi)
         && isA<processorLduInterface>(interfaces[patchi].interface())
        )
        {
            const auto& ppi = refCast<const processorLduInterface>
            (
                interfaces[patchi].interface()
            );

            nbrPerFace[patchi].setSize(localPerFace[patchi].size());
            UIPstream::read
            (
                UPstream::commsTypes::nonBlocking,
                ppi.neighbProcNo(),
                nbrPerFace[patchi],
                ppi.tag(),
                ppi.comm()
            );
        }
    }

    forAll(interfaces, patchi)
    {
        if
        (
            interfaces.set(patchi)
         && isA<processorLduInterface>(interfaces[patchi].interface())
        )
        {
            const auto& ppi = refCast<const processorLduInterface>
            (
                interfaces[patchi].interface()
            );

            UOPstream::write
            (
                UPstream::commsTypes::nonBlocking,
                ppi.neighbProcNo(),
                localPerFace[patchi],
                ppi.tag(),
                ppi.comm()
            );
        }
    }

    UPstream::waitRequests(startOfRequests);

    return nbrPerFace;
}


// Same as exchangePerFaceLabel(), but for a per-interface scalarField
// (fixed-size, one value per face). Used by restrictField()/
// prolongField() to exchange restriction contributions / prolongation
// source values across processor interfaces every V-cycle sweep.
Foam::List<Foam::scalarField> exchangePerFaceScalar
(
    const Foam::lduInterfaceFieldPtrsList& interfaces,
    const Foam::List<Foam::scalarField>& localPerFace
)
{
    using namespace Foam;

    const label startOfRequests = UPstream::nRequests();

    List<scalarField> nbrPerFace(interfaces.size());

    forAll(interfaces, patchi)
    {
        if
        (
            interfaces.set(patchi)
         && isA<processorLduInterface>(interfaces[patchi].interface())
        )
        {
            const auto& ppi = refCast<const processorLduInterface>
            (
                interfaces[patchi].interface()
            );

            nbrPerFace[patchi].setSize(localPerFace[patchi].size());
            UIPstream::read
            (
                UPstream::commsTypes::nonBlocking,
                ppi.neighbProcNo(),
                nbrPerFace[patchi],
                ppi.tag(),
                ppi.comm()
            );
        }
    }

    forAll(interfaces, patchi)
    {
        if
        (
            interfaces.set(patchi)
         && isA<processorLduInterface>(interfaces[patchi].interface())
        )
        {
            const auto& ppi = refCast<const processorLduInterface>
            (
                interfaces[patchi].interface()
            );

            UOPstream::write
            (
                UPstream::commsTypes::nonBlocking,
                ppi.neighbProcNo(),
                localPerFace[patchi],
                ppi.tag(),
                ppi.comm()
            );
        }
    }

    UPstream::waitRequests(startOfRequests);

    return nbrPerFace;
}


// Exchange a per-interface, VARIABLE-length flattened field (e.g. a
// concatenation of variable-size per-face rows), where the receive-side
// total length per interface (nbrTotalSize[patchi]) must already be
// known (typically from a prior exchangePerFaceLabel() row-length
// exchange). Generic over any contiguous Type. Local and remote sizes
// may legitimately DIFFER (that is the whole point of this helper), so
// every receive is posted with its own explicit, correct size BEFORE
// any send is posted - never rely on processorLduInterface::send()'s
// bundled receive-post, which (incorrectly, for this asymmetric case)
// sizes the receive from the LOCAL send buffer.
template<class Type>
Foam::List<Foam::List<Type>> exchangeVariableField
(
    const Foam::lduInterfaceFieldPtrsList& interfaces,
    const Foam::List<Foam::List<Type>>& localFlat,
    const Foam::labelList& nbrTotalSize
)
{
    using namespace Foam;

    const label startOfRequests = UPstream::nRequests();

    List<List<Type>> nbrFlat(interfaces.size());

    forAll(interfaces, patchi)
    {
        if
        (
            interfaces.set(patchi)
         && isA<processorLduInterface>(interfaces[patchi].interface())
        )
        {
            const auto& ppi = refCast<const processorLduInterface>
            (
                interfaces[patchi].interface()
            );

            nbrFlat[patchi].setSize(nbrTotalSize[patchi]);
            UIPstream::read
            (
                UPstream::commsTypes::nonBlocking,
                ppi.neighbProcNo(),
                nbrFlat[patchi],
                ppi.tag(),
                ppi.comm()
            );
        }
    }

    forAll(interfaces, patchi)
    {
        if
        (
            interfaces.set(patchi)
         && isA<processorLduInterface>(interfaces[patchi].interface())
        )
        {
            const auto& ppi = refCast<const processorLduInterface>
            (
                interfaces[patchi].interface()
            );

            UOPstream::write
            (
                UPstream::commsTypes::nonBlocking,
                ppi.neighbProcNo(),
                localFlat[patchi],
                ppi.tag(),
                ppi.comm()
            );
        }
    }

    UPstream::waitRequests(startOfRequests);

    return nbrFlat;
}

} // End anonymous namespace


// * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * //

Foam::List<Foam::DynamicList<Foam::Tuple2<Foam::label, Foam::scalar>>>
Foam::RugeStubenCoarsening::buildRowAdjacency(const lduMatrix& matrix)
{
    const lduAddressing& addr = matrix.lduAddr();
    const label nCells = addr.size();
    const labelUList& owner = addr.lowerAddr();
    const labelUList& neighbour = addr.upperAddr();
    const scalarField& upperCoeffs = matrix.upper();

    List<DynamicList<Tuple2<label, scalar>>> rows(nCells);

    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];
        const scalar coeff = upperCoeffs[facei];

        rows[own].append(Tuple2<label, scalar>(nei, coeff));
        rows[nei].append(Tuple2<label, scalar>(own, coeff));
    }

    return rows;
}


void Foam::RugeStubenCoarsening::calcStrongNeighbours
(
    const lduMatrix& matrix,
    const List<DynamicList<Tuple2<label, scalar>>>& rowAdjacency
)
{
    const label nCells = rowAdjacency.size();
    const scalarField& diagCoeffs = matrix.diag();

    strongNeighbours_.setSize(nCells);

    forAll(rowAdjacency, i)
    {

        // Row measure: largest magnitude of an off-diagonal with the
        // opposite sign to the diagonal. The "textbook" Ruge-Stuben
        // strength test (-a_ij >= theta*max_k(-a_ik)) implicitly assumes
        // a_ii > 0 with (mostly) negative off-diagonals, e.g. as
        // produced by fvm::div/fvm::ddt-dominated matrices. However
        // fvm::laplacian (used for e.g. simpleFoam's pEqn) assembles
        // the mirror-image convention instead: negative diagonal,
        // positive off-diagonals. Testing the off-diagonal sign
        // relative to the row's own diagonal sign - rather than
        // hardcoding "negative" - makes the strength test convention-
        // agnostic and avoids every row spuriously reporting zero
        // strong neighbours (and hence zero coarsening) whenever
        // a_ii < 0.
        const scalar sgnDiag = (diagCoeffs[i] < 0 ? -1 : 1);

        scalar Mi = 0;
        for (const auto& t : rowAdjacency[i])
        {
            const scalar bij = -sgnDiag*t.second();
            if (bij > Mi)
            {
                Mi = bij;
            }
        }

        DynamicList<label> strong;

        if (Mi > 0)
        {
            for (const auto& t : rowAdjacency[i])
            {
                const label j = t.first();
                const scalar bij = -sgnDiag*t.second();

                if (bij >= strengthThreshold_*Mi)
                {
                    strong.append(j);
                }
            }

            Foam::sort(strong);
        }

        strongNeighbours_[i].transfer(strong);
    }
}


void Foam::RugeStubenCoarsening::calcPointType()
{
    const label nCells = strongNeighbours_.size();

    // Transpose graph: influences[j] = { i : j in S_i }
    List<DynamicList<label>> influences(nCells);
    forAll(strongNeighbours_, i)
    {
        for (const label j : strongNeighbours_[i])
        {
            influences[j].append(i);
        }
    }

    // Initial measure lambda_i = |S_i^T|
    labelList lambda(nCells);
    forAll(influences, i)
    {
        lambda[i] = influences[i].size();
    }

    enum class Status : uint8_t { UNDECIDED = 0, IS_C = 1, IS_F = 2 };

    List<Status> status(nCells, Status::UNDECIDED);
    label nUndecided = nCells;

    // ---- Pass 1: greedy measure-based selection ----
    // O(n log n) lazy-invalidated max-priority queue.
    // Entry: (lambda, -cellIndex).  std::priority_queue is a max-heap, so
    // the pair with the largest lambda is always on top; for equal lambdas
    // the smallest cell index wins (matching the original scan-order
    // tie-break because -i is largest when i is smallest).
    // Stale entries (lambda changed after push, or cell already decided)
    // are detected on pop and discarded.
    using PQEntry = std::pair<label, label>; // (lambda, -cellIndex)
    std::priority_queue<PQEntry> pq;
    forAll(lambda, i)
    {
        pq.push({lambda[i], -i});
    }

    while (nUndecided > 0)
    {
        // Pop stale entries until we reach a valid undecided cell.
        label besti = -1;
        while (!pq.empty())
        {
            const auto [lam, negi] = pq.top();
            pq.pop();
            const label i = -negi;
            if (status[i] == Status::UNDECIDED && lambda[i] == lam)
            {
                besti = i;
                break;
            }
        }

        if (besti == -1) break; // queue exhausted – all remaining undecided cells have equal measure 0

        status[besti] = Status::IS_C;
        --nUndecided;

        // F-ify undecided points that strongly depend on besti
        for (const label j : influences[besti])
        {
            if (status[j] == Status::UNDECIDED)
            {
                status[j] = Status::IS_F;
                --nUndecided;

                for (const label k : strongNeighbours_[j])
                {
                    if (status[k] == Status::UNDECIDED)
                    {
                        ++lambda[k];
                        pq.push({lambda[k], -k});
                    }
                }
            }
        }

        // Reduce the measure of undecided points besti depends on
        for (const label k : strongNeighbours_[besti])
        {
            if (status[k] == Status::UNDECIDED)
            {
                --lambda[k];
                pq.push({lambda[k], -k});
            }
        }
    }

    pointType_.setSize(nCells);
    forAll(status, i)
    {
        pointType_[i] =
            (status[i] == Status::IS_C ? pointType::C : pointType::F);
    }

    // ---- Pass 2: F-F common strong-C-neighbour rule ----
    applyFFCommonCNeighbourFixup();

    // ---- Build fine-to-coarse map ----
    fineToCoarse_.setSize(nCells, -1);
    nCoarseCells_ = 0;

    forAll(pointType_, i)
    {
        if (pointType_[i] == pointType::C)
        {
            fineToCoarse_[i] = nCoarseCells_++;
        }
    }
}


void Foam::RugeStubenCoarsening::applyFFCommonCNeighbourFixup()
{
    // For every strongly-connected pair of F-points, ensure they share
    // at least one common strong C-neighbour; otherwise promote one of
    // them to C. Iterated to a fixed point (bounded by nCells passes).
    const label nCells = pointType_.size();

    bool changed = true;
    label safety = 0;

    while (changed && safety++ < nCells)
    {
        changed = false;

        forAll(pointType_, i)
        {
            if (pointType_[i] != pointType::F) continue;

            for (const label j : strongNeighbours_[i])
            {
                if (pointType_[j] != pointType::F) continue;

                bool common = false;
                for (const label c : strongNeighbours_[i])
                {
                    if
                    (
                        pointType_[c] == pointType::C
                     && strongNeighbours_[j].found(c)
                    )
                    {
                        common = true;
                        break;
                    }
                }

                if (!common)
                {
                    pointType_[i] = pointType::C;
                    changed = true;
                    break;
                }
            }
        }
    }
}


void Foam::RugeStubenCoarsening::applyFFCommonCNeighbourFixupParallel
(
    const List<DynamicList<Tuple2<label, label>>>& cellBndFaces,
    const lduInterfaceFieldPtrsList& interfaces,
    const label comm
)
{
    // Bring pointType_ to a purely-local fixed point first - handles
    // every F-F pair whose common C-neighbour (if any) is local.
    applyFFCommonCNeighbourFixup();

    const label nCells = pointType_.size();
    const globalIndex fineGI(nCells, comm);

    // Ghost fine-cell GLOBAL index seen through every boundary face.
    // Purely topological (does not depend on C/F classification), so
    // computed/exchanged once, outside the round loop below.
    List<labelList> localGID(interfaces.size());
    forAll(interfaces, patchi)
    {
        if
        (
            !interfaces.set(patchi)
         || !isA<processorLduInterface>(interfaces[patchi].interface())
        )
        {
            continue;
        }

        const labelUList& faceCells = interfaces[patchi].interface().faceCells();
        localGID[patchi].setSize(faceCells.size());

        forAll(faceCells, facei)
        {
            localGID[patchi][facei] = fineGI.toGlobal(faceCells[facei]);
        }
    }

    const List<labelList> nbrGID(exchangePerFaceLabel(interfaces, localGID));

    // Case 1 grouping: for every ghost GLOBAL index reached via a
    // strong boundary connection, the (possibly several) LOCAL cells
    // that reach it. Two different local cells connected to the SAME
    // ghost (e.g. a ghost bordering more than one of this rank's
    // cells) can serve as each other's cross-rank common C-neighbour
    // with no further communication - purely topological, computed
    // once outside the round loop.
    Map<DynamicList<label>> ghostGIDToLocalCells(2*nCells + 1);
    forAll(cellBndFaces, i)
    {
        for (const auto& pf : cellBndFaces[i])
        {
            const label gid = nbrGID[pf.first()][pf.second()];
            ghostGIDToLocalCells(gid).append(i);
        }
    }

    bool changedGlobal = true;
    label safety = 0;
    const label maxRounds = nCells + 2;

    while (changedGlobal)
    {
        // Refresh ghost C/F status every round: another rank's fixup
        // may have promoted a ghost since the last round.
        List<labelList> localIsC(interfaces.size());
        forAll(interfaces, patchi)
        {
            if
            (
                !interfaces.set(patchi)
             || !isA<processorLduInterface>(interfaces[patchi].interface())
            )
            {
                continue;
            }

            const labelUList& faceCells =
                interfaces[patchi].interface().faceCells();
            localIsC[patchi].setSize(faceCells.size());

            forAll(faceCells, facei)
            {
                localIsC[patchi][facei] =
                (
                    pointType_[faceCells[facei]] == pointType::C ? 1 : 0
                );
            }
        }

        const List<labelList> nbrIsC(exchangePerFaceLabel(interfaces, localIsC));

        // Case 2 exchange: for every boundary face, send the sending
        // rank's own LOCAL strong C-neighbour GLOBAL indices, so the
        // receiving side can check whether any of them is ALSO one of
        // its own (any-interface) ghost connections - i.e. a "shared
        // corner" cell acting as a common C-neighbour for both ranks.
        List<labelList> rowLen(interfaces.size());
        forAll(interfaces, patchi)
        {
            if
            (
                !interfaces.set(patchi)
             || !isA<processorLduInterface>(interfaces[patchi].interface())
            )
            {
                continue;
            }

            const labelUList& faceCells =
                interfaces[patchi].interface().faceCells();
            rowLen[patchi].setSize(faceCells.size());

            forAll(faceCells, facei)
            {
                const label i = faceCells[facei];

                label n = 0;
                if (pointType_[i] == pointType::F)
                {
                    for (const label c : strongNeighbours_[i])
                    {
                        if (pointType_[c] == pointType::C)
                        {
                            ++n;
                        }
                    }
                }
                rowLen[patchi][facei] = n;
            }
        }

        const List<labelList> nbrRowLen(exchangePerFaceLabel(interfaces, rowLen));

        List<labelList> rowStart(interfaces.size());
        labelList nbrTotalSize(interfaces.size(), Zero);
        List<List<label>> localFlat(interfaces.size());

        forAll(interfaces, patchi)
        {
            if
            (
                !interfaces.set(patchi)
             || !isA<processorLduInterface>(interfaces[patchi].interface())
            )
            {
                continue;
            }

            const labelUList& faceCells =
                interfaces[patchi].interface().faceCells();

            rowStart[patchi].setSize(faceCells.size() + 1);
            rowStart[patchi][0] = 0;
            forAll(faceCells, facei)
            {
                rowStart[patchi][facei+1] =
                    rowStart[patchi][facei] + nbrRowLen[patchi][facei];
            }
            nbrTotalSize[patchi] = rowStart[patchi].last();

            label totalLocal = 0;
            forAll(faceCells, facei)
            {
                totalLocal += rowLen[patchi][facei];
            }

            localFlat[patchi].setSize(totalLocal);

            label pos = 0;
            forAll(faceCells, facei)
            {
                const label i = faceCells[facei];
                if (pointType_[i] != pointType::F) continue;

                for (const label c : strongNeighbours_[i])
                {
                    if (pointType_[c] == pointType::C)
                    {
                        localFlat[patchi][pos++] = fineGI.toGlobal(c);
                    }
                }
            }
        }

        const List<List<label>> nbrFlat
        (
            exchangeVariableField<label>(interfaces, localFlat, nbrTotalSize)
        );

        // Build a lookup of every ghost (any interface) this rank can
        // currently see, by GLOBAL index, so Case 2 can be tested with
        // a simple hash lookup rather than a further exchange.
        Map<bool> myGhostIsC(2*nCells + 1);
        forAll(interfaces, patchi)
        {
            if
            (
                !interfaces.set(patchi)
             || !isA<processorLduInterface>(interfaces[patchi].interface())
            )
            {
                continue;
            }

            forAll(nbrGID[patchi], facei)
            {
                myGhostIsC.insert
                (
                    nbrGID[patchi][facei],
                    nbrIsC[patchi][facei] != 0
                );
            }
        }

        bool changedLocal = false;

        forAll(cellBndFaces, i)
        {
            if (pointType_[i] != pointType::F) continue;

            for (const auto& pf : cellBndFaces[i])
            {
                const label patchi = pf.first();
                const label facei = pf.second();

                if (nbrIsC[patchi][facei]) continue; // ghost is C: not F-F

                bool common = false;

                // ---- Case 1 ----
                const label jGID = nbrGID[patchi][facei];
                if (ghostGIDToLocalCells.found(jGID))
                {
                    for (const label c : ghostGIDToLocalCells[jGID])
                    {
                        if (c != i && pointType_[c] == pointType::C)
                        {
                            common = true;
                            break;
                        }
                    }
                }

                // ---- Case 2 ----
                if (!common)
                {
                    const label s = rowStart[patchi][facei];
                    const label e = rowStart[patchi][facei+1];

                    for (label k = s; k < e; ++k)
                    {
                        const label gid = nbrFlat[patchi][k];

                        if
                        (
                            myGhostIsC.found(gid)
                         && myGhostIsC[gid]
                        )
                        {
                            common = true;
                            break;
                        }
                    }
                }

                if (!common)
                {
                    pointType_[i] = pointType::C;
                    changedLocal = true;
                    break;
                }
            }
        }

        // Re-run the local-only pass: promotions above may resolve or
        // create purely-local F-F violations elsewhere.
        if (changedLocal)
        {
            applyFFCommonCNeighbourFixup();
        }

        changedGlobal = returnReduce(changedLocal, orOp<bool>(), comm);

        if (safety++ > maxRounds)
        {
            FatalErrorInFunction
                << "Cross-rank F-F common-C-neighbour fixup failed to"
                   " converge within " << maxRounds << " round(s)."
                << exit(FatalError);
        }
    }
}


void Foam::RugeStubenCoarsening::calcBoundaryStrength
(
    const lduMatrix& matrix,
    const List<DynamicList<Tuple2<label, scalar>>>& rowAdjacency,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces
)
{
    const label nCells = rowAdjacency.size();
    const scalarField& diagCoeffs = matrix.diag();

    // Per-row max magnitude of an opposite-sign-to-diagonal
    // off-diagonal, extended to also consider processor-boundary
    // coupling coefficients (folded in using the same "equivalent
    // internal coefficient" sign convention established in
    // RugeStubenAMG::buildGlobalMatrix(): coeff = -bouCoeffs). This is
    // a local, per-row quantity computed independently of
    // calcStrongNeighbours()'s own (local-only) row measure, and is
    // used only to classify boundary faces here - it does not alter
    // strongNeighbours_.
    scalarField Mi(nCells, Zero);

    forAll(rowAdjacency, i)
    {
        const scalar sgnDiag = (diagCoeffs[i] < 0 ? -1 : 1);

        scalar m = 0;
        for (const auto& t : rowAdjacency[i])
        {
            const scalar bij = -sgnDiag*t.second();
            if (bij > m)
            {
                m = bij;
            }
        }

        Mi[i] = m;
    }

    // A cell may touch more than one processor interface (e.g. at a
    // decomposition corner), so first fold every boundary face's
    // contribution into Mi before testing any of them against the
    // threshold.
    forAll(interfaces, patchi)
    {
        if
        (
            !interfaces.set(patchi)
         || !isA<processorLduInterface>(interfaces[patchi].interface())
        )
        {
            continue;
        }

        const labelUList& faceCells = interfaces[patchi].interface().faceCells();
        const scalarField& bouCoeffs = interfaceBouCoeffs[patchi];

        forAll(faceCells, facei)
        {
            const label i = faceCells[facei];
            const scalar sgnDiag = (diagCoeffs[i] < 0 ? -1 : 1);
            const scalar coeff = -bouCoeffs[facei];
            const scalar bij = -sgnDiag*coeff;

            if (bij > Mi[i])
            {
                Mi[i] = bij;
            }
        }
    }

    ownerStrong_.setSize(interfaces.size());
    boundaryCoeff_.setSize(interfaces.size());
    cellBoundaryFaces_.setSize(nCells);

    forAll(interfaces, patchi)
    {
        if
        (
            !interfaces.set(patchi)
         || !isA<processorLduInterface>(interfaces[patchi].interface())
        )
        {
            continue;
        }

        const labelUList& faceCells = interfaces[patchi].interface().faceCells();
        const scalarField& bouCoeffs = interfaceBouCoeffs[patchi];

        ownerStrong_[patchi].setSize(faceCells.size(), false);
        boundaryCoeff_[patchi].setSize(faceCells.size());

        forAll(faceCells, facei)
        {
            const label i = faceCells[facei];
            const scalar sgnDiag = (diagCoeffs[i] < 0 ? -1 : 1);
            const scalar coeff = -bouCoeffs[facei];
            const scalar bij = -sgnDiag*coeff;

            boundaryCoeff_[patchi][facei] = coeff;

            ownerStrong_[patchi][facei] =
                (Mi[i] > 0 && bij >= strengthThreshold_*Mi[i]);

            cellBoundaryFaces_[i].append(Tuple2<label, label>(patchi, facei));
        }
    }
}


void Foam::RugeStubenCoarsening::calcPointTypeParallel
(
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const label comm
)
{
    const label nCells = strongNeighbours_.size();

    // ---- local S^T (transpose graph), exactly as calcPointType() ----
    List<DynamicList<label>> influences(nCells);
    forAll(strongNeighbours_, i)
    {
        for (const label j : strongNeighbours_[i])
        {
            influences[j].append(i);
        }
    }

    // ---- boundary strength: ownerStrong_[patchi][facei] means the
    //      local faceCell strongly depends on its ghost neighbour ----
    //      (persisted as a member by calcBoundaryStrength(), along
    //      with boundaryCoeff_/cellBoundaryFaces_, for later use by
    //      calcInterpolationParallel()/calcCoarseMatrixParallel())
    {
        const List<DynamicList<Tuple2<label, scalar>>> rowAdjacency
        (
            buildRowAdjacency(matrix)
        );

        calcBoundaryStrength
        (
            matrix,
            rowAdjacency,
            interfaceBouCoeffs,
            interfaces
        );
    }

    // Exchange ownerStrong_ (as 0/1 labels) so every rank also learns,
    // per boundary face, whether its ghost neighbour's row strongly
    // depends on it (necessary since the RS strength test is
    // row-relative and not generally symmetric).
    List<labelList> ownerStrongLabel(interfaces.size());
    forAll(interfaces, patchi)
    {
        if
        (
            interfaces.set(patchi)
         && isA<processorLduInterface>(interfaces[patchi].interface())
        )
        {
            ownerStrongLabel[patchi].setSize(ownerStrong_[patchi].size());
            forAll(ownerStrong_[patchi], facei)
            {
                ownerStrongLabel[patchi][facei] =
                    (ownerStrong_[patchi][facei] ? 1 : 0);
            }
        }
    }

    const List<labelList> nbrStrongLabel
    (
        exchangePerFaceLabel(interfaces, ownerStrongLabel)
    );

    // Per-cell list of (patchi, facei) boundary connections that are
    // strong from either side - the "ghost neighbours" used both for
    // the initial measure and for the CLJP round loop below.
    List<DynamicList<Tuple2<label, label>>> cellBndFaces(nCells);

    forAll(interfaces, patchi)
    {
        if
        (
            !interfaces.set(patchi)
         || !isA<processorLduInterface>(interfaces[patchi].interface())
        )
        {
            continue;
        }

        const labelUList& faceCells = interfaces[patchi].interface().faceCells();

        forAll(faceCells, facei)
        {
            const bool strongEdge =
                ownerStrong_[patchi][facei]
             || (nbrStrongLabel[patchi][facei] != 0);

            if (strongEdge)
            {
                cellBndFaces[faceCells[facei]].append
                (
                    Tuple2<label, label>(patchi, facei)
                );
            }
        }
    }

    // ---- static initial measure lambda_i = |S_i^T|, extended with
    //      ghost points that strongly depend on i across a processor
    //      boundary ----
    labelList lambda(nCells);
    forAll(influences, i)
    {
        lambda[i] = influences[i].size();
    }

    forAll(interfaces, patchi)
    {
        if
        (
            !interfaces.set(patchi)
         || !isA<processorLduInterface>(interfaces[patchi].interface())
        )
        {
            continue;
        }

        const labelUList& faceCells = interfaces[patchi].interface().faceCells();

        forAll(faceCells, facei)
        {
            if (nbrStrongLabel[patchi][facei] != 0)
            {
                ++lambda[faceCells[facei]];
            }
        }
    }

    // ---- static weight w_i = lambda_i (dominant) + globally-unique
    //      tiebreak (subordinate). The global cell index is used as the
    //      tiebreak, which guarantees that NO two cells ever share the
    //      same weight. This is the essential correctness property for
    //      the CLJP Luby-style MIS selection: if two ghost-adjacent
    //      undecided cells have equal weight, NEITHER becomes C (the
    //      >= test on both sides yields isMax=false for both), and the
    //      pair stalls. A hash-based tiebreak (prior implementation)
    //      produced exact cross-rank collisions at the birthday-problem
    //      rate, causing O(n_global) rounds of convergence instead of
    //      O(log n). Using the globally-unique index eliminates all ties.
    const globalIndex globalNumbering(nCells, comm);
    const label nCellsGlobal = globalNumbering.size();
    // Guard against overflow: weight = clampedLambda*(nCellsGlobal+1)+globalI
    // must stay within label range. lambdaClamp is the maximum safe lambda.
    const label lambdaClamp =
        max(label(1), labelMax / (nCellsGlobal + 1));

    labelList weight(nCells);
    forAll(weight, i)
    {
        const label globalI = globalNumbering.toGlobal(i);
        weight[i] = min(lambda[i], lambdaClamp) * (nCellsGlobal + 1) + globalI;
    }

    // Exchange the (static) weight once, per boundary face.
    List<labelList> localWeightPerFace(interfaces.size());
    forAll(interfaces, patchi)
    {
        if
        (
            !interfaces.set(patchi)
         || !isA<processorLduInterface>(interfaces[patchi].interface())
        )
        {
            continue;
        }

        const labelUList& faceCells = interfaces[patchi].interface().faceCells();
        localWeightPerFace[patchi].setSize(faceCells.size());

        forAll(faceCells, facei)
        {
            localWeightPerFace[patchi][facei] = weight[faceCells[facei]];
        }
    }

    const List<labelList> nbrWeightPerFace
    (
        exchangePerFaceLabel(interfaces, localWeightPerFace)
    );

    // ---- CLJP round loop ----
    enum class Status : uint8_t { UNDECIDED = 0, IS_C = 1, IS_F = 2 };

    List<Status> status(nCells, Status::UNDECIDED);
    label nUndecidedLocal = nCells;

    List<labelList> nbrStatusPerFace(interfaces.size());
    forAll(interfaces, patchi)
    {
        if
        (
            interfaces.set(patchi)
         && isA<processorLduInterface>(interfaces[patchi].interface())
        )
        {
            nbrStatusPerFace[patchi].setSize
            (
                interfaces[patchi].interface().faceCells().size(),
                label(Status::UNDECIDED)
            );
        }
    }

    auto exchangeStatus =
    [&]()
    {
        List<labelList> localStatusPerFace(interfaces.size());

        forAll(interfaces, patchi)
        {
            if
            (
                !interfaces.set(patchi)
             || !isA<processorLduInterface>(interfaces[patchi].interface())
            )
            {
                continue;
            }

            const labelUList& faceCells =
                interfaces[patchi].interface().faceCells();

            localStatusPerFace[patchi].setSize(faceCells.size());

            forAll(faceCells, facei)
            {
                localStatusPerFace[patchi][facei] =
                    label(status[faceCells[facei]]);
            }
        }

        nbrStatusPerFace = exchangePerFaceLabel(interfaces, localStatusPerFace);
    };

    // Hard safety cap: each round must decide at least one cell globally, so
    // the number of rounds is bounded by the global cell count, not any single
    // rank's local count. nCellsGlobal is already known from globalNumbering.
    const label maxRounds = nCellsGlobal + 2;
    label round = 0;

    while (true)
    {
        // Determine C-candidates using status/weight fixed as of the
        // start of this round (required so that the local-max
        // argument stays conflict-free across ranks: if two
        // ghost-adjacent points were both to become C, at least one
        // side's comparison must have used the other's up-to-date
        // weight, and a strict ">" test cannot hold both ways at
        // once).
        boolList becomesC(nCells, false);

        forAll(status, i)
        {
            if (status[i] != Status::UNDECIDED) continue;

            bool isMax = true;

            for (const label j : strongNeighbours_[i])
            {
                if
                (
                    status[j] == Status::UNDECIDED
                 && weight[j] >= weight[i]
                )
                {
                    isMax = false;
                    break;
                }
            }

            if (isMax)
            {
                for (const label j : influences[i])
                {
                    if
                    (
                        status[j] == Status::UNDECIDED
                     && weight[j] >= weight[i]
                    )
                    {
                        isMax = false;
                        break;
                    }
                }
            }

            if (isMax)
            {
                for (const auto& pf : cellBndFaces[i])
                {
                    const label patchi = pf.first();
                    const label facei = pf.second();

                    if
                    (
                        nbrStatusPerFace[patchi][facei]
                            == label(Status::UNDECIDED)
                     && nbrWeightPerFace[patchi][facei] >= weight[i]
                    )
                    {
                        isMax = false;
                        break;
                    }
                }
            }

            becomesC[i] = isMax;
        }

        forAll(becomesC, i)
        {
            if (becomesC[i])
            {
                status[i] = Status::IS_C;
                --nUndecidedLocal;
            }
        }

        // Refresh ghost status so the F-ification below (and the next
        // round's local-max test) see this round's new C points.
        exchangeStatus();

        forAll(status, i)
        {
            if (status[i] != Status::UNDECIDED) continue;

            bool makeF = false;

            for (const label j : strongNeighbours_[i])
            {
                if (status[j] == Status::IS_C)
                {
                    makeF = true;
                    break;
                }
            }

            if (!makeF)
            {
                for (const auto& pf : cellBndFaces[i])
                {
                    if
                    (
                        nbrStatusPerFace[pf.first()][pf.second()]
                            == label(Status::IS_C)
                    )
                    {
                        makeF = true;
                        break;
                    }
                }
            }

            if (makeF)
            {
                status[i] = Status::IS_F;
                --nUndecidedLocal;
            }
        }

        // Refresh ghost status again to reflect this round's new F
        // points before the next round (or termination check) uses it.
        exchangeStatus();

        const label nUndecidedGlobal =
            returnReduce(nUndecidedLocal, sumOp<label>(), comm);

        if (nUndecidedGlobal == 0)
        {
            break;
        }

        // Only abort if this rank itself still has undecided cells AND the
        // global round count has exceeded the safety limit. A rank with
        // nUndecidedLocal == 0 is merely participating in collective
        // communication while other ranks finish; it must not abort.
        if (round++ > maxRounds && nUndecidedLocal > 0)
        {
            FatalErrorInFunction
                << "Distributed CLJP-style C/F splitting failed to"
                   " converge within " << maxRounds << " round(s) ("
                << nUndecidedLocal << " undecided cell(s) remain on"
                   " this rank)."
                << exit(FatalError);
        }
    }

    pointType_.setSize(nCells);
    forAll(status, i)
    {
        pointType_[i] =
            (status[i] == Status::IS_C ? pointType::C : pointType::F);
    }


    // ---- Pass 2: F-F common strong-C-neighbour rule, extended to
    //      also cover F-F pairs that span a processor interface (see
    //      applyFFCommonCNeighbourFixupParallel() documentation) ----
    applyFFCommonCNeighbourFixupParallel(cellBndFaces, interfaces, comm);

    // ---- Build fine-to-coarse map ----
    fineToCoarse_.setSize(nCells, -1);
    nCoarseCells_ = 0;

    forAll(pointType_, i)
    {
        if (pointType_[i] == pointType::C)
        {
            fineToCoarse_[i] = nCoarseCells_++;
        }
    }
}


void Foam::RugeStubenCoarsening::calcCoarseGlobalIndex
(
    const lduInterfaceFieldPtrsList& interfaces,
    const label comm
)
{
    coarseGlobalIndexPtr_.reset(new globalIndex(nCoarseCells_, comm));
    const globalIndex& gi = *coarseGlobalIndexPtr_;

    // Per boundary face: this rank's coarse global index if the local
    // faceCell is C, -1 otherwise.
    List<labelList> localCoarseGlobalPerFace(interfaces.size());

    forAll(interfaces, patchi)
    {
        if
        (
            !interfaces.set(patchi)
         || !isA<processorLduInterface>(interfaces[patchi].interface())
        )
        {
            continue;
        }

        const labelUList& faceCells = interfaces[patchi].interface().faceCells();
        localCoarseGlobalPerFace[patchi].setSize(faceCells.size());

        forAll(faceCells, facei)
        {
            const label i = faceCells[facei];

            localCoarseGlobalPerFace[patchi][facei] =
            (
                pointType_[i] == pointType::C
              ? gi.toGlobal(fineToCoarse_[i])
              : -1
            );
        }
    }

    boundaryGhostCoarseIndex_ =
        exchangePerFaceLabel(interfaces, localCoarseGlobalPerFace);
}


void Foam::RugeStubenCoarsening::calcInterpolation
(
    const lduMatrix& matrix,
    const List<DynamicList<Tuple2<label, scalar>>>& rowAdjacency
)
{
    const label nCells = pointType_.size();
    const scalarField& diagCoeffs = matrix.diag();

    interpolationAddr_.setSize(nCells);
    interpolationWeights_.setSize(nCells);

    forAll(pointType_, i)
    {
        if (pointType_[i] == pointType::C)
        {
            interpolationAddr_[i] = labelList(1, fineToCoarse_[i]);
            interpolationWeights_[i] = scalarField(1, scalar(1));
            continue;
        }

        // F-point i: classical (Ruge-Stuben 1987) direct interpolation

        Map<scalar> aRow(2*rowAdjacency[i].size());
        for (const auto& t : rowAdjacency[i])
        {
            aRow.insert(t.first(), t.second());
        }

        scalar Pneg = 0, Ppos = 0, Cneg = 0, Cpos = 0;

        for (const auto& t : rowAdjacency[i])
        {
            const label k = t.first();
            const scalar aik = t.second();

            const bool strongC =
                pointType_[k] == pointType::C
             && strongNeighbours_[i].found(k);

            if (aik < 0)
            {
                Pneg += aik;
                if (strongC)
                {
                    Cneg += aik;
                }
            }
            else if (aik > 0)
            {
                Ppos += aik;
                if (strongC)
                {
                    Cpos += aik;
                }
            }
        }

        // If a sign's off-diagonal mass (Pneg/Ppos) has no strongly
        // connected C-neighbour of the same sign to redistribute onto
        // (Cneg/Cpos == 0), that mass must be lumped into the effective
        // diagonal rather than silently dropped - otherwise conservation
        // is broken and the Galerkin (RAP) coarse matrix can end up with
        // an incorrect (even non-positive) diagonal on subsequent
        // coarsening levels, where off-diagonal signs are typically
        // mixed unlike the original (all-negative-offdiagonal) operator.
        scalar diagNeg = diagCoeffs[i];
        scalar diagPos = diagCoeffs[i];

        // Orphaned mass of a given sign (no strong C-neighbour of that
        // sign to redistribute onto) must be folded into the *other*
        // (surviving) branch's effective diagonal, not its own - the
        // branch it would otherwise modify is only ever read when its
        // own Cneg/Cpos is non-zero, so lumping it into its own
        // diagonal has no effect and silently drops the mass instead.
        if (Cneg == 0)
        {
            diagPos += Pneg;
        }

        if (Cpos == 0)
        {
            diagNeg += Ppos;
        }

        DynamicList<label> addr;
        DynamicList<scalar> w;

        for (const label k : strongNeighbours_[i])
        {
            if (pointType_[k] != pointType::C) continue;

            const scalar aik = aRow[k];

            scalar wik = 0;
            if (aik < 0 && Cneg != 0)
            {
                wik = -(Pneg/Cneg)*(aik/diagNeg);
            }
            else if (aik > 0 && Cpos != 0)
            {
                wik = -(Ppos/Cpos)*(aik/diagPos);
            }

            addr.append(fineToCoarse_[k]);
            w.append(wik);
        }

        interpolationAddr_[i].transfer(addr);
        interpolationWeights_[i].transfer(w);
    }
}


void Foam::RugeStubenCoarsening::calcInterpolationParallel
(
    const lduMatrix& matrix,
    const List<DynamicList<Tuple2<label, scalar>>>& rowAdjacency
)
{
    const label nCells = pointType_.size();
    const scalarField& diagCoeffs = matrix.diag();

    interpolationAddr_.setSize(nCells);
    interpolationWeights_.setSize(nCells);
    interpolationGhostAddr_.setSize(nCells);
    interpolationGhostWeights_.setSize(nCells);

    forAll(pointType_, i)
    {
        if (pointType_[i] == pointType::C)
        {
            interpolationAddr_[i] = labelList(1, fineToCoarse_[i]);
            interpolationWeights_[i] = scalarField(1, scalar(1));
            continue;
        }

        // F-point i: classical (Ruge-Stuben 1987) direct interpolation,
        // extended to also account for strong ghost (across-processor-
        // interface) C-neighbours.

        Map<scalar> aRow(2*rowAdjacency[i].size());
        for (const auto& t : rowAdjacency[i])
        {
            aRow.insert(t.first(), t.second());
        }

        scalar Pneg = 0, Ppos = 0, Cneg = 0, Cpos = 0;

        for (const auto& t : rowAdjacency[i])
        {
            const label k = t.first();
            const scalar aik = t.second();

            const bool strongC =
                pointType_[k] == pointType::C
             && strongNeighbours_[i].found(k);

            if (aik < 0)
            {
                Pneg += aik;
                if (strongC)
                {
                    Cneg += aik;
                }
            }
            else if (aik > 0)
            {
                Ppos += aik;
                if (strongC)
                {
                    Cpos += aik;
                }
            }
        }

        // Boundary (processor-interface) contributions: every boundary
        // connection adds to the row's total off-diagonal mass
        // (Pneg/Ppos), matching the local-adjacency treatment above;
        // only STRONG ghost C-neighbours additionally contribute to
        // Cneg/Cpos (the mass to be redistributed onto strong
        // C-neighbours).
        for (const auto& pf : cellBoundaryFaces_[i])
        {
            const label patchi = pf.first();
            const label facei = pf.second();

            const scalar aik = boundaryCoeff_[patchi][facei];
            const bool strongGhostC =
                ownerStrong_[patchi][facei]
             && boundaryGhostCoarseIndex_[patchi][facei] >= 0;

            if (aik < 0)
            {
                Pneg += aik;
                if (strongGhostC)
                {
                    Cneg += aik;
                }
            }
            else if (aik > 0)
            {
                Ppos += aik;
                if (strongGhostC)
                {
                    Cpos += aik;
                }
            }
        }

        scalar diagNeg = diagCoeffs[i];
        scalar diagPos = diagCoeffs[i];

        // Orphaned mass of a given sign (no strong C-neighbour of that
        // sign to redistribute onto) must be folded into the *other*
        // (surviving) branch's effective diagonal, not its own - the
        // branch it would otherwise modify is only ever read when its
        // own Cneg/Cpos is non-zero, so lumping it into its own
        // diagonal has no effect and silently drops the mass instead.
        if (Cneg == 0)
        {
            diagPos += Pneg;
        }

        if (Cpos == 0)
        {
            diagNeg += Ppos;
        }

        DynamicList<label> addr;
        DynamicList<scalar> w;

        for (const label k : strongNeighbours_[i])
        {
            if (pointType_[k] != pointType::C) continue;

            const scalar aik = aRow[k];

            scalar wik = 0;
            if (aik < 0 && Cneg != 0)
            {
                wik = -(Pneg/Cneg)*(aik/diagNeg);
            }
            else if (aik > 0 && Cpos != 0)
            {
                wik = -(Ppos/Cpos)*(aik/diagPos);
            }

            addr.append(fineToCoarse_[k]);
            w.append(wik);
        }

        interpolationAddr_[i].transfer(addr);
        interpolationWeights_[i].transfer(w);

        // Ghost (across-processor-interface) strong C-neighbour
        // contributions, stored separately (keyed by coarse GLOBAL
        // index) since they are not addressable in this level's local
        // coarse numbering.
        DynamicList<label> gaddr;
        DynamicList<scalar> gw;

        for (const auto& pf : cellBoundaryFaces_[i])
        {
            const label patchi = pf.first();
            const label facei = pf.second();

            if (!ownerStrong_[patchi][facei]) continue;

            const label gk = boundaryGhostCoarseIndex_[patchi][facei];
            if (gk < 0) continue; // ghost is F: no direct P entry here

            const scalar aik = boundaryCoeff_[patchi][facei];

            scalar wik = 0;
            if (aik < 0 && Cneg != 0)
            {
                wik = -(Pneg/Cneg)*(aik/diagNeg);
            }
            else if (aik > 0 && Cpos != 0)
            {
                wik = -(Ppos/Cpos)*(aik/diagPos);
            }

            gaddr.append(gk);
            gw.append(wik);
        }

        interpolationGhostAddr_[i].transfer(gaddr);
        interpolationGhostWeights_[i].transfer(gw);
    }
}


void Foam::RugeStubenCoarsening::calcCoarseMatrix
(
    const lduMatrix& matrix,
    const List<DynamicList<Tuple2<label, scalar>>>& rowAdjacency,
    const label comm
)
{
    const scalarField& diagCoeffs = matrix.diag();

    // Sparse per-row accumulation of the Galerkin (RAP) coarse matrix.
    // Contributions are stored keyed by the *raw* ordered (k, l) index
    // pair (k from the "left" interpolation row, l from the "right"
    // interpolation row) - NOT normalised to (min, max). Normalising to
    // (min, max) would silently sum the (k, l) and (l, k) partial sums
    // together even though they are mathematically identical
    // (P^T A P is symmetric because A is symmetric), which would
    // double-count every off-diagonal entry. Because every adjacency
    // pair (i, j) is visited from both directions (rowAdjacency stores
    // both (i, j) and (j, i), and the diagonal self-term's double loop
    // naturally covers both (k, l) and (l, k)), the raw (k, l) entry is
    // already the complete, correctly-scaled Galerkin coarse-matrix
    // entry - only one triangle needs to be read out at the end.
    std::vector<std::map<label, scalar>> coarseRows(nCoarseCells_);

    auto addTriple =
    [&](const labelList& Pi, const scalarField& Wi,
        const scalar aij,
        const labelList& Pj, const scalarField& Wj)
    {
        forAll(Pi, a)
        {
            const label k = Pi[a];
            const scalar wk = Wi[a];
            if (wk == 0) continue;

            forAll(Pj, b)
            {
                const label l = Pj[b];
                const scalar wl = Wj[b];
                if (wl == 0) continue;

                coarseRows[k][l] += wk*aij*wl;
            }
        }
    };

    forAll(pointType_, i)
    {
        const labelList& Pi = interpolationAddr_[i];
        const scalarField& Wi = interpolationWeights_[i];

        // Diagonal contribution (j == i): full double loop over Pi x Pi
        // is correct here (not a half-loop) since both (k, l) and
        // (l, k) are distinct, independently-valid raw entries.
        addTriple(Pi, Wi, diagCoeffs[i], Pi, Wi);

        // Off-diagonal contributions
        for (const auto& t : rowAdjacency[i])
        {
            const label j = t.first();
            const scalar aij = t.second();

            addTriple(Pi, Wi, aij, interpolationAddr_[j], interpolationWeights_[j]);
        }
    }

    scalarField coarseDiag(nCoarseCells_, Zero);
    DynamicList<label> coarseOwner;
    DynamicList<label> coarseNeighbour;
    DynamicList<scalar> coarseUpper;

    for (label k = 0; k < nCoarseCells_; ++k)
    {
        for (const auto& kv : coarseRows[k])
        {
            const label l = kv.first;
            const scalar val = kv.second;

            if (l == k)
            {
                coarseDiag[k] = val;
            }
            else if (l > k)
            {
                // Only read the upper triangle - the (l, k) mirror
                // entry (computed independently above) carries the
                // same value by construction.
                coarseOwner.append(k);
                coarseNeighbour.append(l);
                coarseUpper.append(val);
            }
        }
    }

    labelList ownerList;
    ownerList.transfer(coarseOwner);

    labelList neighbourList;
    neighbourList.transfer(coarseNeighbour);

    coarseMeshPtr_.reset
    (
        new lduPrimitiveMesh
        (
            nCoarseCells_,
            ownerList,
            neighbourList,
            comm,
            true
        )
    );

    coarseMatrixPtr_.reset(new lduMatrix(*coarseMeshPtr_));

    coarseMatrixPtr_->diag(nCoarseCells_) = coarseDiag;

    // Always allocate upper() (even if empty) - some callers such as
    // LUscalarMatrix::convert() call upper() unconditionally and
    // FatalError if neither upper nor lower is allocated at all.
    coarseMatrixPtr_->upper(coarseUpper.size()) = coarseUpper;
}


void Foam::RugeStubenCoarsening::calcCoarseMatrixParallel
(
    const lduMatrix& matrix,
    const List<DynamicList<Tuple2<label, scalar>>>& rowAdjacency,
    const lduInterfaceFieldPtrsList& interfaces,
    const label comm
)
{
    const scalarField& diagCoeffs = matrix.diag();
    const globalIndex& gi = coarseGlobalIndexPtr_();

    // ------------------------------------------------------------------
    // Step 0: exchange full ghost interpolation rows for boundary F-cells
    // (needed below to form cross-rank RAP triple products i-j-k where j
    // is a ghost F-cell). C-cell ghost rows are NOT exchanged here - they
    // are already fully known locally via boundaryGhostCoarseIndex_
    // (a single entry of weight 1, by construction of interpolation).
    // ------------------------------------------------------------------

    List<labelList> rowLenPerFace(interfaces.size());

    forAll(interfaces, patchi)
    {
        if
        (
            !interfaces.set(patchi)
         || !isA<processorLduInterface>(interfaces[patchi].interface())
        )
        {
            continue;
        }

        const labelUList& faceCells = interfaces[patchi].interface().faceCells();
        rowLenPerFace[patchi].setSize(faceCells.size());

        forAll(faceCells, facei)
        {
            const label i = faceCells[facei];
            rowLenPerFace[patchi][facei] =
            (
                pointType_[i] == pointType::C
              ? 0
              : interpolationAddr_[i].size() + interpolationGhostAddr_[i].size()
            );
        }
    }

    const List<labelList> nbrRowLenPerFace
    (
        exchangePerFaceLabel(interfaces, rowLenPerFace)
    );

    // Prefix-sum offsets into the neighbour's flattened row data, per
    // interface/face, so individual ghost rows can be recovered below.
    List<labelList> nbrRowStart(interfaces.size());
    labelList nbrTotalSize(interfaces.size(), Zero);

    List<List<label>> localFlatAddr(interfaces.size());
    List<List<scalar>> localFlatWeight(interfaces.size());

    forAll(interfaces, patchi)
    {
        if
        (
            !interfaces.set(patchi)
         || !isA<processorLduInterface>(interfaces[patchi].interface())
        )
        {
            continue;
        }

        const labelUList& faceCells = interfaces[patchi].interface().faceCells();

        nbrRowStart[patchi].setSize(faceCells.size() + 1);
        nbrRowStart[patchi][0] = 0;
        forAll(faceCells, facei)
        {
            nbrRowStart[patchi][facei+1] =
                nbrRowStart[patchi][facei] + nbrRowLenPerFace[patchi][facei];
        }
        nbrTotalSize[patchi] = nbrRowStart[patchi].last();

        label totalLocal = 0;
        forAll(faceCells, facei)
        {
            totalLocal += rowLenPerFace[patchi][facei];
        }

        localFlatAddr[patchi].setSize(totalLocal);
        localFlatWeight[patchi].setSize(totalLocal);

        label pos = 0;
        forAll(faceCells, facei)
        {
            const label i = faceCells[facei];
            if (pointType_[i] == pointType::C) continue;

            const labelList& la = interpolationAddr_[i];
            const scalarField& lw = interpolationWeights_[i];
            forAll(la, a)
            {
                localFlatAddr[patchi][pos] = gi.toGlobal(la[a]);
                localFlatWeight[patchi][pos] = lw[a];
                ++pos;
            }

            const labelList& ga = interpolationGhostAddr_[i];
            const scalarField& gw = interpolationGhostWeights_[i];
            forAll(ga, a)
            {
                localFlatAddr[patchi][pos] = ga[a];
                localFlatWeight[patchi][pos] = gw[a];
                ++pos;
            }
        }
    }

    const List<List<label>> nbrFlatAddr
    (
        exchangeVariableField<label>(interfaces, localFlatAddr, nbrTotalSize)
    );
    const List<List<scalar>> nbrFlatWeight
    (
        exchangeVariableField<scalar>(interfaces, localFlatWeight, nbrTotalSize)
    );

    // ------------------------------------------------------------------
    // Step A: local Galerkin (RAP) assembly - identical to the serial
    // calcCoarseMatrix() local-adjacency loop, EXTENDED to also fold in
    // "mixed" local/ghost cross terms.
    //
    // Both the diagonal self-term (i, i) and any internal adjacency
    // edge (i, j) are, in general, RAP triple products over i's FULL
    // interpolation row (interpolationAddr_[i] LOCAL entries PLUS
    // interpolationGhostAddr_[i] GLOBAL/remote entries) crossed with
    // j's FULL row. The original code only summed the LOCAL x LOCAL
    // part (Pi x Pj) - any cell i (or j) that is itself a boundary
    // F-cell with a *mixed* row (some local C-neighbours AND some
    // remote/ghost C-neighbours - the common case for a boundary layer
    // cell, not just the "ghost-only" corner case) silently DROPPED the
    // LOCAL x GHOST / GHOST x LOCAL cross contributions entirely - not
    // merely a convergence-rate approximation, but a genuinely missing
    // piece of the coarse diagonal/coupling that can leave some coarse
    // cell under- (or exactly zero-) diagonal, causing divergence.
    // Unlike the boundary-face terms in Step B (which both ranks
    // compute redundantly and later de-duplicate), these cross terms
    // are computed ENTIRELY from data already local to whichever rank
    // owns cell i (diagCoeffs[i]/aij, Pi and its own ghost weights) -
    // no new communication is needed - but since a DIFFERENT local cell
    // could independently contribute to the same (local, remote) coarse
    // pair, these contributions must be SUMMED (not de-duplicated) both
    // locally here and later across ranks (see diagCrossPerInterface
    // exchange below, which differs from unionCoupling's dedup policy).
    // The GHOST x GHOST corner (both sides remote) is handled
    // separately below (see the Ghost x Ghost block and Step C).
    // ------------------------------------------------------------------

    std::vector<std::map<label, scalar>> coarseRows(nCoarseCells_);

    // Genuine 2-hop "3rd rank" corner couplings queued by addGG() and
    // classifyAndAdd() below (both endpoints/ga,gb are GLOBAL coarse
    // indices), to be resolved by the connectivity-free cross-rank
    // exchange in Step C further down.
    DynamicList<label> pendingGA, pendingGB;
    DynamicList<scalar> pendingVal;

    // Populated by Step C: genuine 2-hop "3rd rank" corner couplings,
    // keyed by the OTHER (remote) rank owning the far side, each entry
    // a canonical GLOBAL (lo, hi) coarse-index pair mapped to its
    // accumulated coefficient. Declared here (rather than local to
    // Step C) so it also remains visible further below when extending
    // the coarse-level interface construction to include brand-new
    // inter-rank interfaces for any remaining (previously unconnected)
    // rank pairs.
    std::map<label, std::map<std::pair<label, label>, scalar>>
        newRemoteCoupling;

    typedef std::map<std::pair<label, label>, scalar> couplingMap;
    List<couplingMap> diagCrossPerInterface(interfaces.size());

    Map<label> procToPatch;
    forAll(interfaces, patchi)
    {
        if
        (
            interfaces.set(patchi)
         && isA<processorLduInterface>(interfaces[patchi].interface())
        )
        {
            const auto& ppi = refCast<const processorLduInterface>
            (
                interfaces[patchi].interface()
            );
            procToPatch.insert(ppi.neighbProcNo(), patchi);
        }
    }

    auto addToCross =
    [&](const label localIdx, const label globalIdx, const scalar contrib)
    {
        const label proci = gi.whichProcID(globalIdx);
        const auto iter = procToPatch.find(proci);
        if (iter == procToPatch.end()) return;

        const label patchi = iter.val();
        diagCrossPerInterface[patchi]
        [
            std::make_pair(localIdx, globalIdx)
        ] += contrib;
    };

    auto addTriple =
    [&](const labelList& Pi, const scalarField& Wi,
        const labelList& Gi, const scalarField& GWi,
        const scalar aij,
        const labelList& Pj, const scalarField& Wj,
        const labelList& Gj, const scalarField& GWj)
    {
        // Local x Local -> coarseRows (as before)
        forAll(Pi, a)
        {
            const label k = Pi[a];
            const scalar wk = Wi[a];
            if (wk == 0) continue;

            forAll(Pj, b)
            {
                const label l = Pj[b];
                const scalar wl = Wj[b];
                if (wl == 0) continue;

                coarseRows[k][l] += wk*aij*wl;
            }

            // Local (k, from i) x Ghost (from j)
            forAll(Gj, b)
            {
                const label lg = Gj[b];
                const scalar wl = GWj[b];
                if (wl == 0) continue;

                addToCross(k, lg, wk*aij*wl);
            }
        }

        // NOTE: no separate "Ghost (from i) x Local (from j)" mirror
        // block here. addTriple() is invoked once per DIRECTION of
        // every internal edge (i, j) - once as addTriple(i-row, aij,
        // j-row) when the outer loop is at i, and again, independently,
        // as addTriple(j-row, aji, i-row) when the outer loop reaches j
        // (rowAdjacency stores both directions, and aij == aji since
        // the matrix is symmetric). The "Local (from i) x Ghost (from
        // j)" block above, run from BOTH of those calls, already
        // produces BOTH "Local-i x Ghost-j" (from the i-call) AND
        // "Local-j x Ghost-i" (from the j-call) exactly once each - a
        // separate mirror block here would add the SAME (localIdx,
        // globalIdx) contribution a second time (for the diagonal
        // self-term, where Pi==Pj and Gi==Gj in a SINGLE call, this
        // was especially direct: both blocks iterated the identical
        // Pi x Gi combinations). Unlike coarseRows (a raw, asymmetric
        // map where both triangles must be filled since only one is
        // ever read back out), diagCrossPerInterface is keyed
        // unambiguously by (localIdx, globalIdx) - it needs each
        // combination added exactly once, not twice.

        // Ghost x Ghost: both sides remote - handled separately below
        // (see the "Ghost x Ghost RAP contributions" block, which
        // queues genuine 3-rank corners for Step C's cross-rank
        // exchange rather than dropping them).
    };

    forAll(pointType_, i)
    {
        const labelList& Pi = interpolationAddr_[i];
        const scalarField& Wi = interpolationWeights_[i];
        const labelList& Gi = interpolationGhostAddr_[i];
        const scalarField& GWi = interpolationGhostWeights_[i];

        addTriple(Pi, Wi, Gi, GWi, diagCoeffs[i], Pi, Wi, Gi, GWi);

        for (const auto& t : rowAdjacency[i])
        {
            const label j = t.first();
            const scalar aij = t.second();

            addTriple
            (
                Pi, Wi, Gi, GWi,
                aij,
                interpolationAddr_[j], interpolationWeights_[j],
                interpolationGhostAddr_[j], interpolationGhostWeights_[j]
            );
        }
    }

    // ------------------------------------------------------------------
    // Ghost x Ghost RAP contributions: the diagonal self-term for a
    // cell i whose interpolation row includes remote C-neighbours
    // (interpolationGhostAddr_[i]), and the internal-edge cross term
    // Gi x Gj for an internal-face pair (i, j) that are both local but
    // each also has remote entries. Step A's addTriple() above only
    // folds in Local x Ghost / Ghost x Local combinations - the pure
    // Ghost x Ghost part is skipped there (both "endpoints" resolve to
    // a coarse cell on ANOTHER rank, so it cannot be added to this
    // rank's own coarseRows/coarseDiag at all). This matters most for
    // boundary F-cells whose ENTIRE row is remote (no local C-neighbour
    // at all) - such a cell's diagonal contribution would otherwise be
    // silently dropped in its entirety, which can leave the owning
    // rank's coarse cell critically under-diagonal (seen in practice
    // at the coarsest levels, where a single missing contribution is a
    // large fraction of that cell's total diagonal).
    //
    // Both ga and gb are, by construction, GLOBAL coarse indices never
    // owned by this (discovering) rank (interpolationGhostAddr_ only
    // ever holds remote strong C-neighbours) - so EVERY Ghost x Ghost
    // triple must be told to whichever rank(s) own ga/gb, with no
    // "compute it locally" fast path at all. Rather than a separate
    // per-interface exchange restricted to cases where an existing
    // fine-level interface happens to connect this rank to the (single
    // or shared) owner, every triple - same-owner or genuine 3-rank
    // corner alike - is queued as a "cell donor" lookup addressed
    // purely by GLOBAL coarse index (pendingGA/pendingGB/pendingVal)
    // and resolved by the single connectivity-free cross-rank exchange
    // in Step C below, which requires no pre-existing interface
    // between this rank and ga's/gb's owner(s) at all.
    // ------------------------------------------------------------------

    {
        auto addGG =
        [&](const label ga, const label gb, const scalar contrib)
        {
            if (contrib == 0) return;

            pendingGA.append(ga);
            pendingGB.append(gb);
            pendingVal.append(contrib);
        };

        forAll(pointType_, i)
        {
            const labelList& Gi = interpolationGhostAddr_[i];
            const scalarField& GWi = interpolationGhostWeights_[i];

            if (Gi.empty()) continue;

            // Diagonal self-term: Gi x Gi (upper triangle, b >= a, to
            // avoid double-processing each unordered (a, b) pair).
            forAll(Gi, a)
            {
                const scalar wa = GWi[a];
                if (wa == 0) continue;

                for (label b = a; b < Gi.size(); ++b)
                {
                    const scalar wb = GWi[b];
                    if (wb == 0) continue;

                    addGG(Gi[a], Gi[b], wa*diagCoeffs[i]*wb);
                }
            }

            // Internal-edge cross term: Gi x Gj. rowAdjacency lists
            // each internal edge from BOTH endpoints, so restrict to
            // i < j to process each edge exactly once.
            for (const auto& t : rowAdjacency[i])
            {
                const label j = t.first();
                if (j <= i) continue;

                const scalar aij = t.second();
                const labelList& Gj = interpolationGhostAddr_[j];
                const scalarField& GWj = interpolationGhostWeights_[j];

                forAll(Gi, a)
                {
                    const scalar wa = GWi[a];
                    if (wa == 0) continue;

                    forAll(Gj, b)
                    {
                        const scalar wb = GWj[b];
                        if (wb == 0) continue;

                        addGG(Gi[a], Gj[b], wa*aij*wb);
                    }
                }
            }
        }
    }

    // ------------------------------------------------------------------
    // Step B: cross-rank RAP contributions from boundary (processor-
    // interface) edges. For each boundary edge (i - ghost j) with fine
    // coefficient aij, fold wk*aij*wl into either this rank's local
    // coarseRows (if the resulting coarse column l is local) or a
    // per-interface coupling map keyed by (local k, remote global l)
    // (if l belongs to this interface's neighbour rank).
    //
    // If l resolves to a THIRD rank (a "2-hop corner" - e.g. j is
    // itself a ghost F-cell whose own ghost C-neighbour lives on a
    // different rank again), the (k, l) triple is queued and routed
    // to l's owning rank via a connectivity-free all-to-all exchange
    // (Step C, below), rather than approximated/dropped.
    // ------------------------------------------------------------------

    typedef std::map<std::pair<label, label>, scalar> couplingMap;
    List<couplingMap> couplingPerInterface(interfaces.size());

    forAll(interfaces, patchi)
    {
        if
        (
            !interfaces.set(patchi)
         || !isA<processorLduInterface>(interfaces[patchi].interface())
        )
        {
            continue;
        }

        const auto& ppi = refCast<const processorLduInterface>
        (
            interfaces[patchi].interface()
        );
        const label neighbProcNo = ppi.neighbProcNo();

        const labelUList& faceCells = interfaces[patchi].interface().faceCells();

        auto classifyAndAdd =
        [&](const label k, const scalar wk, const scalar aij,
            const label l, const scalar wl)
        {
            if (wk == 0 || wl == 0) return;

            const scalar contrib = wk*aij*wl;

            if (gi.isLocal(l))
            {
                coarseRows[k][gi.toLocal(l)] += contrib;
            }
            else if (gi.whichProcID(l) == neighbProcNo)
            {
                couplingPerInterface[patchi][std::make_pair(k, l)] += contrib;
            }
            else
            {
                // Genuine 2-hop "3rd rank" corner: l lives on a THIRD
                // rank not directly reachable via this interface.
                // Queue it for the cross-rank exchange below (Step C)
                // instead of lumping into k's own diagonal - this rank
                // (owning k) already knows both sides of the triple,
                // so only l's remote owner needs telling.
                pendingGA.append(gi.toGlobal(k));
                pendingGB.append(l);
                pendingVal.append(contrib);
            }
        };

        forAll(faceCells, facei)
        {
            const label i = faceCells[facei];
            const scalar aij = boundaryCoeff_[patchi][facei];

            const labelList& Pi = interpolationAddr_[i];
            const scalarField& Wi = interpolationWeights_[i];

            // If i has no LOCAL interpolation entries (i.e. i is an
            // F-cell whose entire row is ghost-only), it cannot anchor
            // any local coarse row k here - its full row is still
            // exchanged above (Step 0) for use by whichever rank sees
            // i as ITS ghost.
            if (Pi.empty()) continue;

            const label gC = boundaryGhostCoarseIndex_[patchi][facei];

            if (gC >= 0)
            {
                // Ghost is a C-point: single-entry row (gC, weight 1).
                forAll(Pi, a)
                {
                    classifyAndAdd(Pi[a], Wi[a], aij, gC, 1.0);
                }
            }
            else
            {
                // Ghost is an F-point: use its exchanged flattened row.
                const label start = nbrRowStart[patchi][facei];
                const label end = nbrRowStart[patchi][facei+1];

                for (label b = start; b < end; ++b)
                {
                    const label l = nbrFlatAddr[patchi][b];
                    const scalar wl = nbrFlatWeight[patchi][b];
                    if (wl == 0) continue;

                    forAll(Pi, a)
                    {
                        classifyAndAdd(Pi[a], Wi[a], aij, l, wl);
                    }
                }
            }
        }
    }

    // Step B covers Pi x Pj for the local entries Pi, but its canonical
    // interface union already supplies off-diagonal transpose terms. Add the
    // missing Gi x Pj complement only when both entries address the same
    // global coarse cell. The neighbour row is either an implicit single C
    // entry or the complete exchanged F row. Both ranks process their own
    // direction, as required for the global (k, k) coefficient.
    forAll(interfaces, patchi)
    {
        if
        (
            !interfaces.set(patchi)
         || !isA<processorLduInterface>(interfaces[patchi].interface())
        )
        {
            continue;
        }

        const labelUList& faceCells = interfaces[patchi].interface().faceCells();

        forAll(faceCells, facei)
        {
            const label i = faceCells[facei];
            const labelList& Gi = interpolationGhostAddr_[i];
            const scalarField& GWi = interpolationGhostWeights_[i];

            if (Gi.empty())
            {
                continue;
            }

            const scalar aij = boundaryCoeff_[patchi][facei];
            const label gC = boundaryGhostCoarseIndex_[patchi][facei];
            const label start = nbrRowStart[patchi][facei];
            const label end = nbrRowStart[patchi][facei+1];

            forAll(Gi, a)
            {
                const label ga = Gi[a];
                const scalar wa = GWi[a];
                if (wa == 0) continue;

                if (gC == ga)
                {
                    pendingGA.append(ga);
                    pendingGB.append(gC);
                    pendingVal.append(wa*aij);
                    continue;
                }

                for (label b = start; b < end; ++b)
                {
                    const label gb = nbrFlatAddr[patchi][b];
                    const scalar wb = nbrFlatWeight[patchi][b];
                    if (wb == 0 || gb != ga) continue;

                    pendingGA.append(ga);
                    pendingGB.append(gb);
                    pendingVal.append(wa*aij*wb);
                }
            }
        }
    }

    // ------------------------------------------------------------------
    // Step C: resolve ALL cross-rank RAP contributions queued above as
    // "cell donor" lookups addressed purely by GLOBAL coarse index
    // (pendingGA/pendingGB/pendingVal - populated by the Ghost x Ghost
    // block above, for BOTH its same-owner and genuine 3-rank-corner
    // cases, and by Step B's third-rank case below). Each triple
    // (ga, gb, contrib) needs telling to whichever of ga's/gb's owning
    // ranks did NOT discover it locally, so every affected rank learns
    // of, and can independently reconstruct, its share of the (ga, gb)
    // coarse-matrix coupling - with no pre-existing fine-level
    // processor interface between this rank and ga's/gb's owner(s)
    // required at all. Uses a connectivity-free, all-to-all exchange
    // (Pstream::exchangeConsensus) that needs no advance knowledge of
    // which ranks talk to which.
    // ------------------------------------------------------------------

    {
        const label myProci = UPstream::myProcNo(comm);

        Map<DynamicList<scalar>> sendBufs;

        forAll(pendingGA, idx)
        {
            const label ga = pendingGA[idx];
            const label gb = pendingGB[idx];
            const scalar val = pendingVal[idx];

            const label procA = gi.whichProcID(ga);
            const label procB = gi.whichProcID(gb);

            const label lo = min(ga, gb);
            const label hi = max(ga, gb);

            if (procA == myProci)
            {
                newRemoteCoupling[procB][std::make_pair(lo, hi)] += val;
            }
            else
            {
                DynamicList<scalar>& buf = sendBufs(procA);
                buf.append(scalar(ga));
                buf.append(scalar(gb));
                buf.append(val);
            }

            // Guard against telling the SAME rank twice when ga and gb
            // are owned by the same remote rank (the Ghost x Ghost
            // block's same-owner case) - the procA branch above (or
            // the corresponding receive-side fold below) already
            // covers that rank exactly once.
            if (procB != procA)
            {
                if (procB == myProci)
                {
                    newRemoteCoupling[procA][std::make_pair(lo, hi)] += val;
                }
                else
                {
                    DynamicList<scalar>& buf = sendBufs(procB);
                    buf.append(scalar(ga));
                    buf.append(scalar(gb));
                    buf.append(val);
                }
            }
        }

        Map<List<scalar>> sendBufsFinal;
        forAllIters(sendBufs, iter)
        {
            // Shrink to logical size before transferring: DynamicList
            // stores its CAPACITY (possibly >= logical size) in the
            // underlying List's own size_ member, so a raw move into a
            // plain List would otherwise pick up trailing, uninitialised
            // capacity slots instead of just the appended entries.
            iter.val().shrink();
            sendBufsFinal.insert(iter.key(), List<scalar>(std::move(iter.val())));
        }

        Map<List<scalar>> recvBufs;
        Pstream::exchangeConsensus<List<scalar>, scalar>
        (
            sendBufsFinal,
            recvBufs,
            UPstream::msgType(),
            comm
        );

        forAllConstIters(recvBufs, iter)
        {
            const List<scalar>& buf = iter.val();
            const label n = buf.size()/3;

            for (label e = 0; e < n; ++e)
            {
                const label ga = label(buf[3*e]);
                const label gb = label(buf[3*e+1]);
                const scalar val = buf[3*e+2];

                const bool gaLocal = gi.isLocal(ga);
                const bool gbLocal = gi.isLocal(gb);

                if (gaLocal && gbLocal)
                {
                    // Ghost x Ghost same-owner case: both endpoints
                    // resolve local to me (I am ga's and gb's sole
                    // owner) - fold directly into coarseRows, exactly
                    // like a local x local RAP contribution.
                    const label la = gi.toLocal(ga);
                    const label lb = gi.toLocal(gb);
                    const label lo = min(la, lb);
                    const label hi = max(la, lb);

                    coarseRows[lo][hi] += val;

                    continue;
                }

                // Genuine corner: exactly one of ga, gb resolves local
                // to me - route the coupling to the OTHER (remote)
                // owner, exactly as before.
                const label remoteProc =
                (
                    gaLocal ? gi.whichProcID(gb) : gi.whichProcID(ga)
                );

                const label lo = min(ga, gb);
                const label hi = max(ga, gb);

                newRemoteCoupling[remoteProc][std::make_pair(lo, hi)] += val;
            }
        }
    }

    // Fold newRemoteCoupling entries whose remote rank already has an
    // EXISTING fine-level interface directly into couplingPerInterface
    // for that interface (converting the canonical GLOBAL (lo, hi) key
    // into the (local k, global l) convention used elsewhere) - these
    // then flow through the existing unionCoupling exchange/dedup
    // logic below exactly like any other Step B contribution. Entries
    // for remote ranks with NO existing interface are left as-is, to
    // be turned into brand-new coarse-level interfaces further below.
    for (const auto& rc : newRemoteCoupling)
    {
        const label remoteProc = rc.first;
        if (!procToPatch.contains(remoteProc)) continue;

        const label patchi = procToPatch[remoteProc];

        for (const auto& kv : rc.second)
        {
            const label lo = kv.first.first;
            const label hi = kv.first.second;
            const scalar val = kv.second;

            const label k =
            (
                gi.isLocal(lo) ? gi.toLocal(lo) : gi.toLocal(hi)
            );
            const label l = gi.isLocal(lo) ? hi : lo;

            couplingPerInterface[patchi][std::make_pair(k, l)] += val;
        }
    }

    // ------------------------------------------------------------------
    // Build the local coarse matrix from coarseRows (identical pattern
    // to the serial calcCoarseMatrix()).
    // ------------------------------------------------------------------

    scalarField coarseDiag(nCoarseCells_, Zero);
    DynamicList<label> coarseOwner;
    DynamicList<label> coarseNeighbour;
    DynamicList<scalar> coarseUpper;

    for (label k = 0; k < nCoarseCells_; ++k)
    {
        for (const auto& kv : coarseRows[k])
        {
            const label l = kv.first;
            const scalar val = kv.second;

            if (l == k)
            {
                coarseDiag[k] += val;
            }
            else if (l > k)
            {
                coarseOwner.append(k);
                coarseNeighbour.append(l);
                coarseUpper.append(val);
            }
        }
    }

    labelList ownerList;
    ownerList.transfer(coarseOwner);

    labelList neighbourList;
    neighbourList.transfer(coarseNeighbour);

    coarseMeshPtr_.reset
    (
        new lduPrimitiveMesh
        (
            nCoarseCells_,
            ownerList,
            neighbourList,
            comm,
            true
        )
    );

    coarseMatrixPtr_.reset(new lduMatrix(*coarseMeshPtr_));

    coarseMatrixPtr_->diag(nCoarseCells_) = coarseDiag;

    coarseMatrixPtr_->upper(coarseUpper.size()) = coarseUpper;

    // ------------------------------------------------------------------
    // Exchange each interface's coupling-pair set (canonicalised to
    // GLOBAL (min, max) coarse-index keys) with the neighbour rank and
    // take the UNION. Each rank's local Step B classification above can
    // independently miss a pair (most notably the dropped 2-hop
    // "corner" term documented above), so without this the two ranks
    // sharing an interface could disagree on whether a coarse-level
    // interface exists at all, or on its face count/order - causing a
    // message-size mismatch (MPI_ERR_TRUNCATE) the first time that
    // coarse level communicates. A set union is commutative, so both
    // ranks compute the IDENTICAL union of (their own local set) and
    // (the set received from the neighbour) - i.e. (localA union
    // localB) on both sides - guaranteeing a symmetric decision and an
    // identical face list/order (std::map's key ordering over
    // std::pair<label,label> is exactly the (min, max) lexicographic
    // order used for canonical sorting).
    // ------------------------------------------------------------------

    typedef std::map<std::pair<label, label>, scalar> keyMap;
    List<keyMap> unionCoupling(interfaces.size());

    {
        List<List<label>> localFlatKeys(interfaces.size());
        List<List<scalar>> localFlatVals(interfaces.size());
        List<labelList> localCountPerFace(interfaces.size());

        forAll(interfaces, patchi)
        {
            if
            (
                !interfaces.set(patchi)
             || !isA<processorLduInterface>(interfaces[patchi].interface())
            )
            {
                continue;
            }

            const label n = label(couplingPerInterface[patchi].size());

            localCountPerFace[patchi].setSize(1, n);
            localFlatKeys[patchi].setSize(2*n);
            localFlatVals[patchi].setSize(n);

            label pos = 0;
            for (const auto& kv : couplingPerInterface[patchi])
            {
                const label Kg = gi.toGlobal(kv.first.first);
                const label Lg = kv.first.second;
                const label minKey = min(Kg, Lg);
                const label maxKey = max(Kg, Lg);

                localFlatKeys[patchi][2*pos] = minKey;
                localFlatKeys[patchi][2*pos+1] = maxKey;
                localFlatVals[patchi][pos] = kv.second;
                ++pos;

                unionCoupling[patchi][std::make_pair(minKey, maxKey)] =
                    kv.second;
            }
        }

        const List<labelList> nbrCountPerFace
        (
            exchangePerFaceLabel(interfaces, localCountPerFace)
        );

        labelList nbrTotalKeySize(interfaces.size(), Zero);
        labelList nbrTotalValSize(interfaces.size(), Zero);

        forAll(interfaces, patchi)
        {
            if
            (
                interfaces.set(patchi)
             && isA<processorLduInterface>(interfaces[patchi].interface())
            )
            {
                nbrTotalValSize[patchi] = nbrCountPerFace[patchi][0];
                nbrTotalKeySize[patchi] = 2*nbrCountPerFace[patchi][0];
            }
        }

        const List<List<label>> nbrFlatKeys
        (
            exchangeVariableField<label>
            (
                interfaces, localFlatKeys, nbrTotalKeySize
            )
        );
        const List<List<scalar>> nbrFlatVals
        (
            exchangeVariableField<scalar>
            (
                interfaces, localFlatVals, nbrTotalValSize
            )
        );

        forAll(interfaces, patchi)
        {
            if
            (
                !interfaces.set(patchi)
             || !isA<processorLduInterface>(interfaces[patchi].interface())
            )
            {
                continue;
            }

            const label n = nbrCountPerFace[patchi][0];

            for (label e = 0; e < n; ++e)
            {
                const label minKey = nbrFlatKeys[patchi][2*e];
                const label maxKey = nbrFlatKeys[patchi][2*e+1];
                const auto key = std::make_pair(minKey, maxKey);

                if
                (
                    unionCoupling[patchi].find(key)
                 == unionCoupling[patchi].end()
                )
                {
                    unionCoupling[patchi][key] = nbrFlatVals[patchi][e];
                }
            }
        }
    }

    // ------------------------------------------------------------------
    // Exchange and SUM (not de-duplicate) the diagonal local x ghost
    // cross-term contributions (diagCrossPerInterface, built in Step A
    // above). Unlike unionCoupling's face-based terms, both ranks do
    // NOT compute the same value for a given key here - each rank's
    // contribution comes from its OWN cells' diagonals and is
    // independent of (additive to) whatever the neighbour rank
    // contributes for the same coarse pair - so matching keys must be
    // ADDED together, not de-duplicated.
    // ------------------------------------------------------------------

    List<keyMap> diagCrossFinal(interfaces.size());

    {
        List<List<label>> localFlatKeys(interfaces.size());
        List<List<scalar>> localFlatVals(interfaces.size());
        List<labelList> localCountPerFace(interfaces.size());

        forAll(interfaces, patchi)
        {
            if
            (
                !interfaces.set(patchi)
             || !isA<processorLduInterface>(interfaces[patchi].interface())
            )
            {
                continue;
            }

            const label n = label(diagCrossPerInterface[patchi].size());

            localCountPerFace[patchi].setSize(1, n);
            localFlatKeys[patchi].setSize(2*n);
            localFlatVals[patchi].setSize(n);

            label pos = 0;
            for (const auto& kv : diagCrossPerInterface[patchi])
            {
                const label Kg = gi.toGlobal(kv.first.first);
                const label Lg = kv.first.second;
                const label minKey = min(Kg, Lg);
                const label maxKey = max(Kg, Lg);

                localFlatKeys[patchi][2*pos] = minKey;
                localFlatKeys[patchi][2*pos+1] = maxKey;
                localFlatVals[patchi][pos] = kv.second;
                ++pos;

                diagCrossFinal[patchi][std::make_pair(minKey, maxKey)] +=
                    kv.second;
            }
        }

        const List<labelList> nbrCountPerFace
        (
            exchangePerFaceLabel(interfaces, localCountPerFace)
        );

        labelList nbrTotalKeySize(interfaces.size(), Zero);
        labelList nbrTotalValSize(interfaces.size(), Zero);

        forAll(interfaces, patchi)
        {
            if
            (
                interfaces.set(patchi)
             && isA<processorLduInterface>(interfaces[patchi].interface())
            )
            {
                nbrTotalValSize[patchi] = nbrCountPerFace[patchi][0];
                nbrTotalKeySize[patchi] = 2*nbrCountPerFace[patchi][0];
            }
        }

        const List<List<label>> nbrFlatKeys
        (
            exchangeVariableField<label>
            (
                interfaces, localFlatKeys, nbrTotalKeySize
            )
        );
        const List<List<scalar>> nbrFlatVals
        (
            exchangeVariableField<scalar>
            (
                interfaces, localFlatVals, nbrTotalValSize
            )
        );

        forAll(interfaces, patchi)
        {
            if
            (
                !interfaces.set(patchi)
             || !isA<processorLduInterface>(interfaces[patchi].interface())
            )
            {
                continue;
            }

            const label n = nbrCountPerFace[patchi][0];

            for (label e = 0; e < n; ++e)
            {
                const label minKey = nbrFlatKeys[patchi][2*e];
                const label maxKey = nbrFlatKeys[patchi][2*e+1];
                const auto key = std::make_pair(minKey, maxKey);

                diagCrossFinal[patchi][key] += nbrFlatVals[patchi][e];
            }
        }
    }

    // Merge the (de-duplicated) face-based unionCoupling with the
    // (summed) diagonal cross terms into one final coupling map per
    // interface.
    List<keyMap> combinedCoupling(interfaces.size());
    forAll(unionCoupling, patchi)
    {
        combinedCoupling[patchi] = unionCoupling[patchi];
        for (const auto& kv : diagCrossFinal[patchi])
        {
            combinedCoupling[patchi][kv.first] += kv.second;
        }
    }

    // ------------------------------------------------------------------
    // TEMPORARY DIAGNOSTIC (RSDEBUG-GALERKINCHECK): gather the raw
    // ingredients of the Galerkin product (fine matrix A + interpolation
    // operator P, exactly as fed into Step A/B/GG above) from every rank
    // to the master, independently compute the EXACT reference P^T A P
    // there via a single-rank triple loop with no ghost/interface
    // complexity at all (the same algorithm as the serial
    // calcCoarseMatrix()'s addTriple()), and diff it entry-by-entry
    // against the ACTUAL distributed result (coarseRows +
    // combinedCoupling, gathered the same way). This isolates the
    // correctness of the distributed RAP ASSEMBLY logic from the (also
    // suspect, but separately-varying) C/F splitting algorithm - both
    // sides use the SAME, parallel-derived P. Runs ONCE (first call
    // only, via a static guard) since it is O(nFineGlobal) memory/time
    // on the master and only needed for initial validation.
    // ------------------------------------------------------------------
    {
        static bool validated = false;
        if (!validated)
        {
            validated = true;

            const label nCellsLocal = pointType_.size();
            const globalIndex fineGI(nCellsLocal, comm);
            const label nProcs = UPstream::nProcs(comm);
            const label myProci = UPstream::myProcNo(comm);

            // ---- Fine matrix diagonal ----
            List<scalarField> allDiag(nProcs);
            allDiag[myProci] = diagCoeffs;
            Pstream::gatherList(allDiag, UPstream::msgType(), comm);

            // ---- Fine matrix internal (local-local) off-diagonal edges ----
            DynamicList<label> locI, locJ;
            DynamicList<scalar> locVal;
            forAll(rowAdjacency, i)
            {
                for (const auto& t : rowAdjacency[i])
                {
                    locI.append(fineGI.toGlobal(i));
                    locJ.append(fineGI.toGlobal(t.first()));
                    locVal.append(t.second());
                }
            }
            List<labelList> allEdgeI(nProcs), allEdgeJ(nProcs);
            List<scalarField> allEdgeVal(nProcs);
            allEdgeI[myProci] = locI;
            allEdgeJ[myProci] = locJ;
            allEdgeVal[myProci] = locVal;
            Pstream::gatherList(allEdgeI, UPstream::msgType(), comm);
            Pstream::gatherList(allEdgeJ, UPstream::msgType(), comm);
            Pstream::gatherList(allEdgeVal, UPstream::msgType(), comm);

            // ---- Fine matrix boundary (processor-interface) edges: need
            //      the GHOST's GLOBAL fine index, obtained by exchanging
            //      this rank's own fineGI-converted faceCells across
            //      every processor interface. ----
            List<labelList> localFineGID(interfaces.size());
            forAll(interfaces, patchi)
            {
                if
                (
                    !interfaces.set(patchi)
                 || !isA<processorLduInterface>(interfaces[patchi].interface())
                )
                {
                    continue;
                }

                const labelUList& faceCells =
                    interfaces[patchi].interface().faceCells();
                localFineGID[patchi].setSize(faceCells.size());
                forAll(faceCells, facei)
                {
                    localFineGID[patchi][facei] =
                        fineGI.toGlobal(faceCells[facei]);
                }
            }
            const List<labelList> nbrFineGID
            (
                exchangePerFaceLabel(interfaces, localFineGID)
            );

            DynamicList<label> bndI, bndJ;
            DynamicList<scalar> bndVal;
            forAll(interfaces, patchi)
            {
                if
                (
                    !interfaces.set(patchi)
                 || !isA<processorLduInterface>(interfaces[patchi].interface())
                )
                {
                    continue;
                }

                const labelUList& faceCells =
                    interfaces[patchi].interface().faceCells();
                forAll(faceCells, facei)
                {
                    bndI.append(localFineGID[patchi][facei]);
                    bndJ.append(nbrFineGID[patchi][facei]);
                    bndVal.append(boundaryCoeff_[patchi][facei]);
                }
            }
            List<labelList> allBndI(nProcs), allBndJ(nProcs);
            List<scalarField> allBndVal(nProcs);
            allBndI[myProci] = bndI;
            allBndJ[myProci] = bndJ;
            allBndVal[myProci] = bndVal;
            Pstream::gatherList(allBndI, UPstream::msgType(), comm);
            Pstream::gatherList(allBndJ, UPstream::msgType(), comm);
            Pstream::gatherList(allBndVal, UPstream::msgType(), comm);

            // ---- Interpolation operator P: flatten local + ghost rows,
            //      keyed by GLOBAL fine row and GLOBAL coarse column ----
            DynamicList<label> pRow, pCol;
            DynamicList<scalar> pVal;
            forAll(pointType_, i)
            {
                const label gI = fineGI.toGlobal(i);

                const labelList& Pi = interpolationAddr_[i];
                const scalarField& Wi = interpolationWeights_[i];
                forAll(Pi, a)
                {
                    pRow.append(gI);
                    pCol.append(gi.toGlobal(Pi[a]));
                    pVal.append(Wi[a]);
                }

                const labelList& Gi = interpolationGhostAddr_[i];
                const scalarField& GWi = interpolationGhostWeights_[i];
                forAll(Gi, a)
                {
                    pRow.append(gI);
                    pCol.append(Gi[a]); // already global
                    pVal.append(GWi[a]);
                }
            }
            List<labelList> allPRow(nProcs), allPCol(nProcs);
            List<scalarField> allPVal(nProcs);
            allPRow[myProci] = pRow;
            allPCol[myProci] = pCol;
            allPVal[myProci] = pVal;
            Pstream::gatherList(allPRow, UPstream::msgType(), comm);
            Pstream::gatherList(allPCol, UPstream::msgType(), comm);
            Pstream::gatherList(allPVal, UPstream::msgType(), comm);

            // ---- This rank's own contribution to the ACTUAL distributed
            //      result: coarseRows (local-local, unique per global k)
            //      plus combinedCoupling (per-interface, GLOBAL (min,max)
            //      keyed - gathered from BOTH sides so rank-to-rank
            //      (a)symmetry can be cross-checked on the master) ----
            DynamicList<label> actK, actL;
            DynamicList<scalar> actVal;
            for (label k = 0; k < nCoarseCells_; ++k)
            {
                for (const auto& kv : coarseRows[k])
                {
                    actK.append(gi.toGlobal(k));
                    actL.append(gi.toGlobal(kv.first));
                    actVal.append(kv.second);
                }
            }
            forAll(combinedCoupling, patchi)
            {
                for (const auto& kv : combinedCoupling[patchi])
                {
                    actK.append(kv.first.first);
                    actL.append(kv.first.second);
                    actVal.append(kv.second);
                }
            }
            // Brand-new (Step C) interfaces for remote ranks with no
            // existing fine-level interface are NOT part of
            // combinedCoupling - include them here too, otherwise this
            // diagnostic would blindly report the fix as ineffective.
            for (const auto& rc : newRemoteCoupling)
            {
                if (procToPatch.contains(rc.first)) continue;

                for (const auto& kv : rc.second)
                {
                    actK.append(kv.first.first);
                    actL.append(kv.first.second);
                    actVal.append(kv.second);
                }
            }
            List<labelList> allActK(nProcs), allActL(nProcs);
            List<scalarField> allActVal(nProcs);
            allActK[myProci] = actK;
            allActL[myProci] = actL;
            allActVal[myProci] = actVal;
            Pstream::gatherList(allActK, UPstream::msgType(), comm);
            Pstream::gatherList(allActL, UPstream::msgType(), comm);
            Pstream::gatherList(allActVal, UPstream::msgType(), comm);

            if (UPstream::master(comm))
            {
                // ---- Build global fine A (symmetric, canonical
                //      (min,max) off-diagonal map + diagonal array) ----
                label nFineGlobal = 0;
                for (label p = 0; p < nProcs; ++p)
                {
                    nFineGlobal += allDiag[p].size();
                }

                scalarField ADiag(nFineGlobal, Zero);
                {
                    label offset = 0;
                    for (label p = 0; p < nProcs; ++p)
                    {
                        forAll(allDiag[p], i)
                        {
                            ADiag[offset + i] = allDiag[p][i];
                        }
                        offset += allDiag[p].size();
                    }
                }

                std::map<std::pair<label, label>, scalar> Aoff;
                for (label p = 0; p < nProcs; ++p)
                {
                    forAll(allEdgeI[p], e)
                    {
                        const label gi1 = allEdgeI[p][e];
                        const label gj1 = allEdgeJ[p][e];
                        const label lo = min(gi1, gj1);
                        const label hi = max(gi1, gj1);
                        Aoff[std::make_pair(lo, hi)] = allEdgeVal[p][e];
                    }
                    forAll(allBndI[p], e)
                    {
                        const label gi1 = allBndI[p][e];
                        const label gj1 = allBndJ[p][e];
                        const label lo = min(gi1, gj1);
                        const label hi = max(gi1, gj1);
                        Aoff[std::make_pair(lo, hi)] = allBndVal[p][e];
                    }
                }

                // ---- Build global P: fine global row -> list of
                //      (coarse global col, weight) ----
                std::vector<std::vector<std::pair<label, scalar>>> Prows
                (
                    nFineGlobal
                );
                for (label p = 0; p < nProcs; ++p)
                {
                    forAll(allPRow[p], e)
                    {
                        Prows[allPRow[p][e]].push_back
                        (
                            std::make_pair(allPCol[p][e], allPVal[p][e])
                        );
                    }
                }

                // ---- Reference Galerkin product: C_ref = P^T A P,
                //      computed exactly as the serial calcCoarseMatrix()
                //      addTriple() lambda, but over the FULL global
                //      problem (no ghost/interface distinction at all). ----
                std::map<std::pair<label, label>, scalar> Cref;
                auto addRef =
                [&](const label k, const scalar wk,
                    const scalar aij,
                    const label l, const scalar wl,
                    const char* tag = "")
                {
                    if (wk == 0 || wl == 0) return;
                    Cref[std::make_pair(k, l)] += wk*aij*wl;
                };

                for (label i = 0; i < nFineGlobal; ++i)
                {
                    // diagonal self-term
                    for (const auto& pk : Prows[i])
                    {
                        for (const auto& pl : Prows[i])
                        {
                            addRef
                            (
                                pk.first, pk.second,
                                ADiag[i],
                                pl.first, pl.second,
                                "diag-self"
                            );
                        }
                    }
                }
                for (const auto& kv : Aoff)
                {
                    const label i = kv.first.first;
                    const label j = kv.first.second;
                    const scalar aij = kv.second;

                    for (const auto& pk : Prows[i])
                    {
                        for (const auto& pl : Prows[j])
                        {
                            addRef(pk.first, pk.second, aij, pl.first, pl.second, "edge-fwd");
                        }
                    }
                    for (const auto& pk : Prows[j])
                    {
                        for (const auto& pl : Prows[i])
                        {
                            addRef(pk.first, pk.second, aij, pl.first, pl.second, "edge-mirror");
                        }
                    }
                }


                // Fold Cref to canonical (min,max) keys - it is exactly
                // symmetric by construction (A symmetric), so each
                // off-diagonal physical pair was independently filled
                // TWICE above (once as (k,l), once as (l,k), both with
                // the identical value) - halve after folding. Diagonal
                // (k,k) entries were only ever filled once (the pure
                // diagonal self-term double loop) - leave as-is.
                std::map<std::pair<label, label>, scalar> CrefCanon;
                for (const auto& kv : Cref)
                {
                    const label lo = min(kv.first.first, kv.first.second);
                    const label hi = max(kv.first.first, kv.first.second);
                    CrefCanon[std::make_pair(lo, hi)] += kv.second;
                }
                for (auto& kv : CrefCanon)
                {
                    if (kv.first.first != kv.first.second)
                    {
                        kv.second *= 0.5;
                    }
                }

                // ---- Actual distributed result, gathered the same way.
                //      Insert-else-compare (not +=): every physical
                //      (k,l) pair should arrive with the SAME value from
                //      every direction/rank it is reported from (local
                //      mirror direction, or the neighbour rank sharing
                //      an interface) - any genuine mismatch here is
                //      itself a direct proof of a rank-asymmetric bug. ----
                std::map<std::pair<label, label>, scalar> Cact;
                label nMismatchKey = 0;
                for (label p = 0; p < nProcs; ++p)
                {
                    forAll(allActK[p], e)
                    {
                        const label k = allActK[p][e];
                        const label l = allActL[p][e];
                        const label lo = min(k, l);
                        const label hi = max(k, l);
                        const auto key = std::make_pair(lo, hi);

                        const auto it = Cact.find(key);
                        if (it == Cact.end())
                        {
                            Cact[key] = allActVal[p][e];
                        }
                        else if (mag(it->second - allActVal[p][e]) > 1e-9)
                        {
                            ++nMismatchKey;
                        }
                    }
                }

                Pout<< "RSDEBUG-GALERKINCHECK nFineGlobal=" << nFineGlobal
                    << " nCoarseGlobal=" << gi.size()
                    << " nRankSymmetryMismatches=" << nMismatchKey
                    << endl;

                scalar maxAbsDiff = 0;
                scalar maxRelDiff = 0;
                label nCompared = 0;
                label nDiffOver1pct = 0;
                std::pair<label, label> worstKey(-1, -1);

                std::map<std::pair<label, label>, scalar> allKeys;
                for (const auto& kv : CrefCanon) allKeys[kv.first] = 0;
                for (const auto& kv : Cact) allKeys[kv.first] = 0;

                for (const auto& kv : allKeys)
                {
                    const scalar refVal =
                    (
                        CrefCanon.count(kv.first)
                      ? CrefCanon[kv.first] : scalar(0)
                    );
                    const scalar actVal =
                    (
                        Cact.count(kv.first)
                      ? Cact[kv.first] : scalar(0)
                    );

                    const scalar diff = mag(actVal - refVal);
                    const scalar scale = max(mag(refVal), SMALL);
                    const scalar rel = diff/scale;

                    ++nCompared;
                    if (diff > maxAbsDiff) maxAbsDiff = diff;
                    if (rel > maxRelDiff)
                    {
                        maxRelDiff = rel;
                        worstKey = kv.first;
                    }
                    if (rel > 0.01) ++nDiffOver1pct;
                }

                Pout<< "RSDEBUG-GALERKINCHECK nCompared=" << nCompared
                    << " nDiffOver1pct=" << nDiffOver1pct
                    << " maxAbsDiff=" << maxAbsDiff
                    << " maxRelDiff=" << maxRelDiff
                    << " worstKey=(" << worstKey.first << ","
                    << worstKey.second << ")"
                    << " refVal="
                    << (CrefCanon.count(worstKey) ? CrefCanon[worstKey] : 0)
                    << " actVal="
                    << (Cact.count(worstKey) ? Cact[worstKey] : 0)
                    << endl;

                // Dump every mismatching entry, along with which global
                // rank each coarse cell k/l belongs to (via gi), to look
                // for a structural pattern (e.g. "3-rank chain" pairs
                // that are missed entirely by the local-only combinedCoupling
                // exchange logic).
                for (const auto& kv : allKeys)
                {
                    const scalar refVal =
                    (
                        CrefCanon.count(kv.first)
                      ? CrefCanon[kv.first] : scalar(0)
                    );
                    const scalar actVal =
                    (
                        Cact.count(kv.first)
                      ? Cact[kv.first] : scalar(0)
                    );
                    const scalar diff = mag(actVal - refVal);
                    const scalar scale = max(mag(refVal), SMALL);
                    if (diff/scale > 0.01)
                    {
                        Pout<< "RSDEBUG-GALERKINCHECK-MISS k="
                            << kv.first.first
                            << " (rank " << gi.whichProcID(kv.first.first)
                            << ") l=" << kv.first.second
                            << " (rank " << gi.whichProcID(kv.first.second)
                            << ") refVal=" << refVal
                            << " actVal=" << actVal
                            << endl;
                    }
                }
            }
        }
    }

    DynamicList<label> validPatchIs;
    forAll(combinedCoupling, patchi)
    {
        if (!combinedCoupling[patchi].empty())
        {
            validPatchIs.append(patchi);
        }
    }

    // Genuine 2-hop "3rd rank" corner couplings (Step C, above) whose
    // remote rank does NOT already have an existing fine-level
    // interface need a brand-new coarse-level interface, built below
    // alongside (but separately from) the validPatchIs-derived ones.
    DynamicList<label> newRemoteProcs;
    for (const auto& rc : newRemoteCoupling)
    {
        if (!procToPatch.contains(rc.first))
        {
            newRemoteProcs.append(rc.first);
        }
    }

    const label nTotalIfaces = validPatchIs.size() + newRemoteProcs.size();

    coarseInterfaceFields_.setSize(nTotalIfaces);
    coarseInterfaces_.setSize(nTotalIfaces);
    coarseInterfaceBouCoeffs_.setSize(nTotalIfaces);

    // The new processorGAMGInterface objects built below are handed to
    // coarseMeshPtr_->addInterfaces() further down, which makes the
    // mesh's own (owning) primitiveInterfaces_ the SOLE owner (see the
    // comment below). Collect them here in a local, non-owning
    // lduInterfacePtrsList - NOT in a separately-owned PtrList member -
    // to avoid a double-free (previously: coarseGAMGInterfaces_ AND
    // coarseMeshPtr_'s primitiveInterfaces_ both believed they owned
    // the same objects, crashing with a use-after-free/double-delete
    // in ~lduPrimitiveMesh() at RugeStubenCoarsening destruction time).
    lduInterfacePtrsList tempInterfaces(nTotalIfaces);

    forAll(validPatchIs, vi)
    {
        const label patchi = validPatchIs[vi];

        const auto& ppi = refCast<const processorLduInterface>
        (
            interfaces[patchi].interface()
        );

        // combinedCoupling is keyed by canonical (min, max) GLOBAL
        // coarse index pairs, and std::map iterates keys in ascending
        // (lexicographic) order - i.e. already exactly the desired
        // (min, max) canonical face order, identically on both ranks.
        // No additional sort is needed.
        std::vector<std::pair<std::pair<label, label>, scalar>> entries
        (
            combinedCoupling[patchi].begin(),
            combinedCoupling[patchi].end()
        );

        labelList faceCells(entries.size());
        scalarField coeffs(entries.size());
        forAll(entries, fi)
        {
            const label minKey = entries[fi].first.first;
            const label maxKey = entries[fi].first.second;

            // Exactly one of {minKey, maxKey} is a coarse GLOBAL index
            // owned by this rank (the other belongs to neighbProcNo) -
            // recover this rank's LOCAL coarse cell index for the new
            // interface's faceCells, regardless of whether this key
            // was discovered locally or only via the neighbour.
            faceCells[fi] =
            (
                gi.isLocal(minKey)
              ? gi.toLocal(minKey)
              : gi.toLocal(maxKey)
            );

            // entries[fi].second ("contrib") was accumulated using
            // boundaryCoeff_ (the NEGATED raw bouCoeffs, i.e. the
            // "internal-equivalent" sign convention - see
            // calcPointTypeParallel()'s "coeff = -bouCoeffs[facei]").
            // That convention is correct for coarseDiag/coarseUpper
            // (genuine internal matrix entries, built from the same
            // contrib just above), but coarseInterfaceBouCoeffs_ feeds
            // GaussSeidelSmoother/lduMatrix::residual()/Amul() via
            // processorGAMGInterfaceField::updateInterfaceMatrix(),
            // which expect the RAW bouCoeffs convention (the opposite
            // sign) - negate here to convert back.
            coeffs[fi] = -entries[fi].second;
        }

        tempInterfaces.set
        (
            vi,
            new processorGAMGInterface
            (
                vi,
                lduInterfacePtrsList(0),
                faceCells,
                identity(faceCells.size()),
                comm,
                ppi.myProcNo(),
                ppi.neighbProcNo(),
                tensorField(0),
                ppi.tag()
            )
        );

        coarseInterfaceFields_.set
        (
            vi,
            new processorGAMGInterfaceField
            (
                refCast<const processorGAMGInterface>(tempInterfaces[vi]),
                false,
                0
            )
        );

        coarseInterfaces_.set(vi, &coarseInterfaceFields_[vi]);
        coarseInterfaceBouCoeffs_.set(vi, new scalarField(std::move(coeffs)));
    }

    // Brand-new interfaces for genuine 2-hop "3rd rank" corner
    // couplings (Step C, above) whose remote rank has no existing
    // fine-level interface. Both ranks of such a pair always compute
    // an IDENTICAL entry set (Step C's exchange delivers the same
    // (ga, gb, contrib) triples to BOTH owning ranks) and, for the
    // message tag, a deterministic formula from (this rank, remote
    // rank) that both sides evaluate identically - no negotiation
    // round needed. Each interface is self-contained and identified
    // only by its own myProcNo/neighbProcNo/tag/comm (exactly as for
    // the validPatchIs-derived interfaces above), so the two ranks do
    // not need to agree on a common patch/vi ordering.
    {
        const label myProci = UPstream::myProcNo(comm);
        const label nProcsGlobal = UPstream::nProcs(comm);

        forAll(newRemoteProcs, ni)
        {
            const label vi = validPatchIs.size() + ni;
            const label remoteProc = newRemoteProcs[ni];

            std::vector<std::pair<std::pair<label, label>, scalar>> entries
            (
                newRemoteCoupling[remoteProc].begin(),
                newRemoteCoupling[remoteProc].end()
            );

            labelList faceCells(entries.size());
            scalarField coeffs(entries.size());
            forAll(entries, fi)
            {
                const label minKey = entries[fi].first.first;
                const label maxKey = entries[fi].first.second;

                faceCells[fi] =
                (
                    gi.isLocal(minKey)
                  ? gi.toLocal(minKey)
                  : gi.toLocal(maxKey)
                );

                // Same sign-convention negation as validPatchIs above.
                coeffs[fi] = -entries[fi].second;
            }

            const label pairTag =
                1000000
              + min(myProci, remoteProc)*nProcsGlobal
              + max(myProci, remoteProc);

            tempInterfaces.set
            (
                vi,
                new processorGAMGInterface
                (
                    vi,
                    lduInterfacePtrsList(0),
                    faceCells,
                    identity(faceCells.size()),
                    comm,
                    myProci,
                    remoteProc,
                    tensorField(0),
                    pairTag
                )
            );

            coarseInterfaceFields_.set
            (
                vi,
                new processorGAMGInterfaceField
                (
                    refCast<const processorGAMGInterface>(tempInterfaces[vi]),
                    false,
                    0
                )
            );

            coarseInterfaces_.set(vi, &coarseInterfaceFields_[vi]);
            coarseInterfaceBouCoeffs_.set(vi, new scalarField(std::move(coeffs)));
        }
    }

    // The coarse lduPrimitiveMesh above was built via the "components,
    // without interfaces" constructor (owner/neighbour/comm only) - its
    // OWN lduAddressing therefore has ZERO patches registered until we
    // explicitly add them here. Without this, lduAddr().patchAddr() is
    // an empty list, and any code that indexes it per the (non-empty)
    // coarseInterfaces_ built above (e.g. GaussSeidelSmoother::smooth(),
    // lduMatrix::residual()/Amul()) dereferences out-of-range and
    // FatalErrors ("Cannot dereference nullptr at index 0 in range
    // [0,0)"). lduPrimitiveMesh::addInterfaces() takes logical
    // ownership of the objects referenced by tempInterfaces (copies the
    // pointers into its own OWNING primitiveInterfaces_ PtrList) - see
    // lduPrimitiveMesh.C. tempInterfaces itself (a non-owning UPtrList)
    // does NOT delete anything when it goes out of scope, so ownership
    // cleanly transfers exactly once, with coarseMeshPtr_ as the sole
    // owner for the remaining lifetime of this object.
    coarseMeshPtr_->addInterfaces
    (
        tempInterfaces,
        lduPrimitiveMesh::nonBlockingSchedule<processorGAMGInterface>
        (
            tempInterfaces
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RugeStubenCoarsening::RugeStubenCoarsening
(
    const lduMatrix& matrix,
    const scalar strengthThreshold,
    const label comm
)
:
    strengthThreshold_(strengthThreshold),
    pointType_(),
    strongNeighbours_(),
    fineToCoarse_(),
    nCoarseCells_(0),
    interpolationAddr_(),
    interpolationWeights_(),
    coarseMeshPtr_(),
    coarseMatrixPtr_()
{
    if (!matrix.symmetric())
    {
        FatalErrorInFunction
            << "RugeStubenCoarsening currently only supports symmetric"
               " (SPD) matrices"
            << exit(FatalError);
    }

    const List<DynamicList<Tuple2<label, scalar>>> rowAdjacency
    (
        buildRowAdjacency(matrix)
    );

    calcStrongNeighbours(matrix, rowAdjacency);
    calcPointType();
    calcInterpolation(matrix, rowAdjacency);
    calcCoarseMatrix(matrix, rowAdjacency, comm);
}


Foam::RugeStubenCoarsening::RugeStubenCoarsening
(
    const lduMatrix& matrix,
    const scalar strengthThreshold,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const label comm
)
:
    strengthThreshold_(strengthThreshold),
    pointType_(),
    strongNeighbours_(),
    fineToCoarse_(),
    nCoarseCells_(0),
    interpolationAddr_(),
    interpolationWeights_(),
    coarseMeshPtr_(),
    coarseMatrixPtr_()
{
    if (!matrix.symmetric())
    {
        FatalErrorInFunction
            << "RugeStubenCoarsening currently only supports symmetric"
               " (SPD) matrices"
            << exit(FatalError);
    }

    const List<DynamicList<Tuple2<label, scalar>>> rowAdjacency
    (
        buildRowAdjacency(matrix)
    );

    calcStrongNeighbours(matrix, rowAdjacency);
    calcPointTypeParallel(matrix, interfaceBouCoeffs, interfaces, comm);
    calcCoarseGlobalIndex(interfaces, comm);
    calcInterpolationParallel(matrix, rowAdjacency);
    calcCoarseMatrixParallel(matrix, rowAdjacency, interfaces, comm);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::RugeStubenCoarsening::restrictField
(
    scalarField& coarseFld,
    const scalarField& fineFld,
    const lduInterfaceFieldPtrsList& fineInterfaces
) const
{
    coarseFld.setSize(nCoarseCells_);
    coarseFld = Zero;

    // Local contributions (identical to the purely-serial case)
    forAll(pointType_, i)
    {
        const labelList& addr = interpolationAddr_[i];
        const scalarField& w = interpolationWeights_[i];

        forAll(addr, a)
        {
            coarseFld[addr[a]] += w[a]*fineFld[i];
        }
    }

    if (cellBoundaryFaces_.empty())
    {
        // Plain (serial) constructor was used - no boundary bookkeeping
        return;
    }

    // Cross-rank ghost contributions: for every boundary face whose
    // local cell strongly depends on a ghost C-neighbour, send this
    // cell's contribution (interpolation weight * fineFld[i]) across
    // the interface; 0 for faces with no such contribution. Fixed-size
    // per-face exchange - both ranks of every interface always send
    // and receive, so this is symmetric/deadlock-free.
    List<scalarField> sendVal(fineInterfaces.size());

    forAll(fineInterfaces, patchi)
    {
        if
        (
            !fineInterfaces.set(patchi)
         || !isA<processorLduInterface>(fineInterfaces[patchi].interface())
        )
        {
            continue;
        }

        const labelUList& faceCells =
            fineInterfaces[patchi].interface().faceCells();

        sendVal[patchi].setSize(faceCells.size(), Zero);
    }

    forAll(pointType_, i)
    {
        if (pointType_[i] == pointType::C) continue;

        label pos = 0;
        for (const auto& pf : cellBoundaryFaces_[i])
        {
            const label patchi = pf.first();
            const label facei = pf.second();

            if
            (
                !ownerStrong_[patchi][facei]
             || boundaryGhostCoarseIndex_[patchi][facei] < 0
            )
            {
                continue;
            }

            sendVal[patchi][facei] =
                interpolationGhostWeights_[i][pos]*fineFld[i];

            ++pos;
        }
    }

    const List<scalarField> recvVal
    (
        exchangePerFaceScalar(fineInterfaces, sendVal)
    );

    forAll(fineInterfaces, patchi)
    {
        if
        (
            !fineInterfaces.set(patchi)
         || !isA<processorLduInterface>(fineInterfaces[patchi].interface())
        )
        {
            continue;
        }

        const labelUList& faceCells =
            fineInterfaces[patchi].interface().faceCells();

        forAll(faceCells, facei)
        {
            const label j = faceCells[facei];

            if (fineToCoarse_[j] >= 0)
            {
                coarseFld[fineToCoarse_[j]] += recvVal[patchi][facei];
            }
        }
    }
}


void Foam::RugeStubenCoarsening::prolongField
(
    scalarField& fineFld,
    const scalarField& coarseFld,
    const lduInterfaceFieldPtrsList& fineInterfaces
) const
{
    fineFld.setSize(pointType_.size());

    List<scalarField> sendVal(fineInterfaces.size());

    if (!cellBoundaryFaces_.empty())
    {
        // Every boundary face sends its own local cell's coarse value
        // (if it is a C-point), so the neighbouring rank can use it as
        // a ghost prolongation source; 0 (unused) otherwise.
        forAll(fineInterfaces, patchi)
        {
            if
            (
                !fineInterfaces.set(patchi)
             || !isA<processorLduInterface>(fineInterfaces[patchi].interface())
            )
            {
                continue;
            }

            const labelUList& faceCells =
                fineInterfaces[patchi].interface().faceCells();

            sendVal[patchi].setSize(faceCells.size(), Zero);

            forAll(faceCells, facei)
            {
                const label j = faceCells[facei];

                if (pointType_[j] == pointType::C)
                {
                    sendVal[patchi][facei] = coarseFld[fineToCoarse_[j]];
                }
            }
        }
    }

    const List<scalarField> recvVal
    (
        cellBoundaryFaces_.empty()
      ? List<scalarField>(fineInterfaces.size())
      : exchangePerFaceScalar(fineInterfaces, sendVal)
    );

    forAll(pointType_, i)
    {
        const labelList& addr = interpolationAddr_[i];
        const scalarField& w = interpolationWeights_[i];

        scalar val = 0;
        forAll(addr, a)
        {
            val += w[a]*coarseFld[addr[a]];
        }

        if (!cellBoundaryFaces_.empty() && pointType_[i] != pointType::C)
        {
            label pos = 0;
            for (const auto& pf : cellBoundaryFaces_[i])
            {
                const label patchi = pf.first();
                const label facei = pf.second();

                if
                (
                    !ownerStrong_[patchi][facei]
                 || boundaryGhostCoarseIndex_[patchi][facei] < 0
                )
                {
                    continue;
                }

                val +=
                    interpolationGhostWeights_[i][pos]
                   *recvVal[patchi][facei];

                ++pos;
            }
        }

        fineFld[i] = val;
    }
}


// ************************************************************************* //
