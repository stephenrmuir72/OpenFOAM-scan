/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2024 OpenCFD Ltd.
    Copyright (C) 2026 Keysight Technologies
    Copyright (C) 2026 Mark Olesen
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

Class
    Foam::distributedTriSurfaceMesh

Description
    IOobject and searching on distributed triSurface. All processors hold
    (possibly overlapping) part of the overall surface. All queries are
    distributed to the processor that can answer it and the result sent back.

    Can work in different modes:
    - frozen : no changes or updating
    - follow : makes sure each processor has all the triangles inside the
    externally provided bounding box (usually the mesh bounding box).
    Guarantees minimum amount of communication since mesh-local queries
    should be answerable without any comms.
    - independent : surface is decomposed according to the triangle centres
    so the decomposition might be radically different from the mesh
    decomposition. Guarantees best memory balance but at the expense of
    more communication.
    - distributed : has already been distributed and can be used as is.
    - nodeMaster : the node-local master has the complete surface and answers
      queries. No communication between nodes.
    - fullTree : this is like \em distributed/independent
      with the surface kd-tree decomposed equally onto each rank.
    - sparseTree : this is like \em distributed/independent
      with the surface kd-tree decomposed across node masters

    Note: addressing used:
    distributedTriSurfaceMesh: none

    triSurfaceMesh:
    - surf.pointFaces()     : edge addressing (for volume tests only)
    - surf.edges()          : edgeTree
    - surf.faceFaces()      : only if minQuality > 0

Usage
    In e.g. \c system/snappyHexMeshDict:
    \verbatim
    {
        myGeometry
        {
            type distributedTriSurfaceMesh;
            file "box.obj";

            numberOfSubdomains  8;
            method              hierarchical;
            n                   (2 2 2);
        }
    }
    \endverbatim

    \verbatim
        // Mandatory (inherited) entries
        type            distributedTriSurfaceMesh;

        // Optional entries

        //- distributionType (default is distributed)
        distributionType        nodeMaster;

        //- Use bounding boxes (default) or unique decomposition of triangles
        //- (i.e. do not duplicate triangles)
        decomposeUsingBbs       false;

        //- For nodeMaster/sparseTree mode : override the number of ranks
        //- per master (default: all ranks on the same node)
        nProcessorsPerMaster    4;

        //- Calculate bb that does not overlap any of the other bounding boxes.
        exclusionBb             false;

        //- Merging distance for points on the surface. Only used when reading.
        mergeDistance           1e-10;

        //- Override the closed state of the initial surface. Only used when
        //  reading.
        closed                  false;

        //- Override the outside volume type. Only used when reading.
        outsideVolumeType       outside;
    }
    \endverbatim

SourceFiles
    distributedTriSurfaceMesh.cxx
    distributedTriSurfaceMesh.txx
    distributedTriSurfaceMeshGeometry.cxx

\*---------------------------------------------------------------------------*/

#ifndef Foam_distributedTriSurfaceMesh_H
#define Foam_distributedTriSurfaceMesh_H

#include "triSurfaceMesh.H"
#include "localIOdictionary.H"
#include "IOdictionary.H"
#include "Pair.H"
#include "globalIndex.H"
#include "DynamicField.H"
#include "boundBoxTree.H"
#include "triangleFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class mapDistribute;
class decompositionMethod;

// Typedefs
typedef Pair<point> segment;


/*---------------------------------------------------------------------------*\
                   Class distributedTriSurfaceMesh Declaration
\*---------------------------------------------------------------------------*/

class distributedTriSurfaceMesh
:
    public triSurfaceMesh
{
public:

    // Static Data

        //- The type of distribution style to use
        enum distributionType : int8_t
        {
            FROZEN = -1,   //!< placeholder type
            FOLLOW = 0,
            TRANSPARENT,
            INDEPENDENT,
            DISTRIBUTED,
            NODEMASTER,
            FULLTREE,      //!< surface parts on each rank
            SPARSETREE,    //!< surface parts on each node master
        };

        static const Enum<distributionType> distributionTypeNames_;

        //- Order in which receives surfaces are processed by
        //- internal routines (eg, exchangeSurfaces)
        enum class receiveOrder
        {
            sequenced,   //!< In rank order. Receives are still all pre-posted
            unsequenced  //!< As they complete. Must be order-independent
        };

private:

    // Data Types

        //- Enumerated control or (advanced) tuning parameters
        //- [normally used as part of a bit mask]
        enum controlType : unsigned
        {
            //- Use direct send/recv instead of PstreamBuffers for
            //- surface exchange
            ctrlDirectExchange = 1,

            //- Overlap communication and querying in node-master mode
            ctrlOverlappingNodeMaster = 2,

            //- Use polling (instead of PstreamBuffers) for
            //- handling surface overlap calculations
            ctrlPollingOverlaps = 4,

            //- Default preferences (compile-time)
            ctrlDefault =
            (
                ctrlDirectExchange
              | ctrlOverlappingNodeMaster
              | ctrlPollingOverlaps
            )
        };


    // Private Member Data

        mutable autoPtr<IOdictionary> decomposeParDict_;

        //- Decomposition used when independently decomposing surface.
        mutable autoPtr<decompositionMethod> decomposer_;

        //- Bounding box settings
        localIOdictionary dict_;

        //- Bounding boxes of all processors, with tree searching
        boundBoxTree procBb_;

        //- if !decomposeUsingBbs_ (so no fill in but unique processor
        //  assignment) we can still have a local boundbox where we guarantee
        //  we have the interacting triangles. This is an optional filtering
        //  at start up. Is list over all processors - currently only
        //  a single bb per processor.
        boundBoxTree exclusionBb_;

        //- Global triangle numbering
        mutable autoPtr<globalIndex> globalTris_;

        //- Optional per-vertex normals. TBD: move to triSurface? or have
        //- per-triangle 3 normals so we can interpolate and have features
        mutable autoPtr<List<FixedList<vector, 3>>> vertexNormals_;

        //- Use bounding boxes (default) or unique decomposition of triangles
        //- (i.e. do not duplicate triangles)
        bool decomposeUsingBbs_;

        //- Calculate bb that does not overlap any of the other bounding
        //- boxes
        bool doExclusionBb_;

        //- The (wanted) distribution type.
        distributionType distType_;

        //- The (current) distribution type. Used to trigger re-distribution
        //  when starting from undecomposed surface.
        distributionType currentDistType_;

        //- Merging distance
        scalar mergeDist_;

        //- Control/tuning parameters
        unsigned control_;

        //- Optional communicator to send queries. Used in NODEMASTER mode.
        int comm_;

        //- In nodeMaster mode : IDs of master processors
        //- (one per node, or as specified by nProcessorsPerMaster)
        DynamicList<int> masterProcIDs_;


    // Private Member Functions

        //- Running in "transparent mode" (passes straight to triSurfaceMesh)
        //- * for non-parallel
        //- * if flagged internally as "transparent"
        //- * running NODEMASTER, but with 1 rank per master (COMM_SELF)
        //- .
        bool runningTransparent() const noexcept
        {
            return
            (
                !UPstream::parRun() ||
                distType_ == distributionType::TRANSPARENT ||
                (
                    distType_ == distributionType::NODEMASTER
                 && comm_ == UPstream::commSelf()
                )
            );
        }

        //- Frees the comm_ communicator and reset to -1.
        //- Does not attempt to free any persistent communicators
        void freeCommunicator()
        {
            UPstream::freeLocalCommunicator(comm_);
            comm_ = -1;
        }

        //- Adapt to changes in dictionary syntax between versions
        static void upgradeDictionary(dictionary& dict);


        // Read

            //- Search for io.local directory (=triSurface) either in case
            //  directory or in parent directory
            static word findLocalInstance(const IOobject& io);

            //- Read my additional data from dictionary. Additional flag to
            //- say whether we can use master-only geometric tests.
            bool readSettings(bool isUndecomposed);

            //- [FUTURE] Adjust decomposition settings according to
            //- a given heuristic.
            bool optimizeStrategy();

            //- Define/redefine the masterProcIDs_ and adjust the
            //- the communicator (comm_)
            void setNodeMasters
            (
                int nProcessorsPerMaster,
                bool needsCommunicator
            );

            //- Construction helper: generate vertex normals upon reading
            //- undecomposed surface
            void calcVertexNormals
            (
                const triSurface& surf,
                List<FixedList<vector, 3>>& vn
            ) const;


        // Line intersection

            //- 'subtract' bounding box \p b from \p a such that
            //- max volume remains
            static treeBoundBox subtract
            (
                const treeBoundBox& a,
                const treeBoundBox& b
            );

            //- Update exclusionBb_ from new procBb.
            void calcExclusionBb();

            //- Is segment local to myProcNo
            bool isInExclusion(const point& start, const point& end) const;

            //- Split segment into subsegments overlapping the processor
            //  bounding box.
            //void Foam::distributedTriSurfaceMesh::splitSegment
            //(
            //    const label segmentI,
            //    const point& start,
            //    const point& end,
            //    const treeBoundBox& bb,
            //
            //    DynamicList<segment>& allSegments,
            //    DynamicList<label>& allSegmentMap,
            //    DynamicList<label> sendMap
            //) const

            //- Distribute segments into overlapping processor
            //  bounding boxes. Sort per processor.
            void distributeSegment
            (
                const label,
                const point& start,
                const point& end,

                DynamicList<segment>&,
                DynamicList<label>&,
                List<DynamicList<label>>&
            ) const;

            //- Divide edges into local and remote segments. Construct map to
            //  distribute and collect data.
            autoPtr<mapDistribute> distributeSegments
            (
                const pointField& start,
                const pointField& end,
                const bitSet& doTest,

                List<segment>& allSegments,
                List<label>& allSegmentMap
            ) const;

            //- Split edges, distribute, test and collect.
            void findLine
            (
                const bool nearestIntersection,
                const pointField& start,
                const pointField& end,
                List<pointIndexHit>& info
            ) const;


        // Triangle index

            //- Helper: inplace convert local triangle indices to global ones
            void convertTriIndices
            (
                UList<pointIndexHit>& info
            ) const;

            //- Obtains global indices from pointIndexHit and swaps them back
            //- to their original processor.
            //- Used to calculate local region and normal.
            autoPtr<mapDistribute> calcLocalQueries
            (
                const UList<pointIndexHit>&,
                labelList& triangleIndex
            ) const;


        // Nearest

            //- Calculate map to send centres+radius to processors
            autoPtr<mapDistribute> calcLocalQueries
            (
                const bool includeLocalProcessor,
                const pointField& centres,
                const scalarField& radiusSqr,

                pointField& allCentres,
                scalarField& allRadiusSqr,
                labelList& allSegmentMap
            ) const;


        // Side

            //- Side of nearest point on faces w.r.t. samples.
            //- Handles nearest on edge/point
            void surfaceSide
            (
                const pointField& samples,
                const UList<pointIndexHit>& nearestInfo,
                List<volumeType>& region
            ) const;


            // Caching of volume type (based on indexedOctree)

            //- Set node type on any node containing the triangle
            volumeType markMixed
            (
                const indexedOctree<treeDataTriSurface>& tree,
                const label nodei,
                const triPointRef& tri,
                PackedList<2>& nodeTypes
            ) const;

            //- Set node type on any node containing triangles from
            //- the surface
            void markMixed
            (
                const indexedOctree<treeDataTriSurface>& tree,
                const label nodei,
                const triSurface& surf,
                PackedList<2>& nodeTypes
            ) const;

            //- Set node type on any node overlapping any remote triangles.
            //- Only valid if using unique decomposition.
            void markMixedOverlap
            (
                const indexedOctree<treeDataTriSurface>& tree,
                PackedList<2>& nodeTypes
            ) const;

            //- Collect mid points of tree boxes
            void collectLeafMids
            (
                const label nodeI,
                DynamicField<point>& midPoints
            ) const;

            //- Find volume type of tree boxes
            volumeType calcVolumeType
            (
                const UList<volumeType>& midPointTypes,
                label& index,
                PackedList<2>& nodeTypes,
                const label nodeI
            ) const;

            //- Calculate inside/outside of midpoint of tree nodes
            void cacheVolumeType(PackedList<2>& nt) const;

            //- Look up any cached data. Return unknown if cannot be determined.
            volumeType cachedVolumeType
            (
                const label nodeI,
                const point& sample
            ) const;


        // Surface redistribution

            //- Helper: get decompositionMethod
            const decompositionMethod& decomposer() const;

            //- Finds new bounds based on an independent decomposition.
            void calcDistribution
            (
                const triSurface& s,
                //! [out] the decomposition (depending on distributionType)
                labelList& distribution
            ) const;

            //- Find triangle \p f in allFaces. Returns -1 if not found
            static label findTriangle
            (
                //! The faces to search (haystack)
                const UList<labelledTri>& allFaces,

                //! The point to face connectivity information
                const UList<labelList>& allPointFaces,

                //! The face to find (needle)
                const labelledTri& f
            );

            //- Find triangle \p f in allFaces. Returns -1 if not found
            static label findTriangle
            (
                //! The faces to search (haystack)
                const UList<labelledTri>& allFaces,

                //! The point to face connectivity information
                //! indexed by lowest vertex point
                const Map<DynamicList<label>>& allPointFaces,

                //! The face to find (needle)
                const labelledTri& f
            );

            //- Merge triSurface (subTris, subPoints) into allTris, allPoints.
            static void merge
            (
                const scalar mergeDist,
                const UList<labelledTri>& subTris,
                const UList<point>& subPoints,

                DynamicList<labelledTri>& allTris,
                DynamicList<point>& allPoints,

                labelList& faceConstructMap,
                labelList& pointConstructMap
            );

            //- Distribute stored fields
            template<class Type>
            void distributeFields(const mapDistribute& map);

            //- Calculate the distribution of triangles based on the
            //- distribution type.
            //  Gets passed in mesh bounding boxes (if in 'follow' mode)
            //
            //  \returns true if the \p distribution is in use.
            bool getDistribution
            (
                const distributionType distType,
                // The bounding box to follow (for FOLLOW mode)
                const UList<treeBoundBox>& followBb,
                //! [out] the decomposition (depending on distributionType)
                labelList& distribution,
                //! [out] the proc bounding box of all (destination) triangles
                List<treeBoundBox>& newProcBb
            ) const;

            //- Do actual sending and receiving of triangles and points.
            void distributeGeometry
            (
                labelListList&& faceSendMap,
                autoPtr<mapDistribute>& faceMap,
                autoPtr<mapDistribute>& pointMap,
                List<labelledTri>& allTris,
                pointField& allPoints
            ) const;

            //- Do actual sending and receiving of triangles and points.
            //- Using PstreamBuffers for the exchange
            void distributeGeometry_buffered
            (
                labelListList&& faceSendMap,
                autoPtr<mapDistribute>& faceMap,
                autoPtr<mapDistribute>& pointMap,
                List<labelledTri>& allTris,
                pointField& allPoints
            ) const;


        // Low-level communication

            //- What type of data exchange is needed?
            enum class sendStyles : int
            {
                mode_general,   //!< multiple ranks send
                mode_decompose, //!< only master sends
                mode_broadcast, //!< master sends identical content
                mode_identity   //!< special case of broadcast (no-subset)
            };

            //- Determine which type of exchange is needed based on
            //- the face send maps
            static sendStyles analyzeSendMaps
            (
                const triSurface& surf,
                const UList<labelList>& faceSendMap,
                int communicator
            );

            //- Send surface primitives (default: non-blocking)
            static void sendSurfaceGeometry
            (
                int toProcNo,           //!< Send to this rank
                const UList<point>& points,
                const UList<labelledTri>& faces,
                int baseMessageTag,     //!< The base message tag
                int communicator,
                UPstream::commsTypes ctype = UPstream::commsTypes::nonBlocking
            );

            //- Send surface primitives (default: non-blocking)
            static void sendSurfaceGeometry
            (
                int toProcNo,           //!< Send to this rank
                const triSurface& s,
                int baseMessageTag,     //!< The base message tag
                int communicator,
                UPstream::commsTypes ctype = UPstream::commsTypes::nonBlocking
            )
            {
                sendSurfaceGeometry
                (
                    toProcNo,
                    s.points(),
                    s.surfFaces(),
                    baseMessageTag,
                    communicator,
                    ctype
                );
            }

            //- Receive surface primitives (default: non-blocking)
            static void recvSurfaceGeometry
            (
                int fromProcNo,         //!< Receive from this rank
                List<point>& points,
                List<labelledTri>& faces,
                int baseMessageTag,     //!< The base message tag
                int communicator,
                UPstream::commsTypes ctype = UPstream::commsTypes::nonBlocking
            );

            //- Receive surface primitives (default: non-blocking)
            static void recvSurfaceGeometry
            (
                int fromProcNo,         //!< Receive from this rank
                triSurface& s,
                int baseMessageTag,     //!< The base message tag
                int communicator,
                UPstream::commsTypes ctype = UPstream::commsTypes::nonBlocking
            )
            {
                recvSurfaceGeometry
                (
                    fromProcNo,
                    s.storedPoints(),
                    s.storedFaces(),
                    baseMessageTag,
                    communicator,
                    ctype
                );
            }

            //- Sparse point-to-point surface exchange using the default
            //- sendSurfaceGeometry/recvSurfaceGeometry (non-blocking).
            //  Wrapper for the more granular exchangeSurfaces helper
            //  (see distributedTriSurfaceMesh.C).
            //  - sendOp(proci)  -> autoPtr<triSurface> (null: nothing sent)
            //  - procOp(proci, triSurface&) consumes a completed receive
            //  - preRecvOp(const UList<labelPair>&) optional, called with the
            //    per-rank receive sizes (nPoints, nFaces) before receiving
            template
            <
                class SendOp,
                class ProcessOp,
                class PreRecvOp = Foam::identityOp
            >
            static void exchangeSurfaces
            (
                const int communicator,
                const receiveOrder recvOrder,
                const int baseMessageTag,
                SendOp&& sendOp,
                ProcessOp&& procOp,
                PreRecvOp&& preRecvOp = PreRecvOp()
            );

public:

    // Generated Methods

        //- No copy construct
        distributedTriSurfaceMesh(const distributedTriSurfaceMesh&) = delete;

        //- No copy assignment
        void operator=(const distributedTriSurfaceMesh&) = delete;


    //- Runtime type information
    TypeName("distributedTriSurfaceMesh");


    // Constructors

        //- Construct from components. Assumes triSurface is already decomposed
        //- and dictionary contains corresponding information
        distributedTriSurfaceMesh
        (
            const IOobject&,
            const triSurface&,
            const dictionary& dict
        );

        //- Construct read. Does findInstance to find io.local()
        //-   - if found local : assume distributed
        //-   - if found in parent : assume undistributed. Can e.g. check for
        //-     closedness.
        distributedTriSurfaceMesh(const IOobject& io);

        //- Construct from dictionary (used by searchableSurface).
        //  Does read. Does findInstance to find io.local().
        //-   - if found local : assume distributed
        //-   - if found in parent : assume undistributed. Can e.g. check for
        //-     closedness.
        distributedTriSurfaceMesh
        (
            const IOobject& io,
            const dictionary& dict
        );


    //- Destructor
    virtual ~distributedTriSurfaceMesh();

        //- Clear storage
        void clearOut();


    // Member Functions

        //- Triangle indexing (demand driven)
        const globalIndex& globalTris() const;


        // Methods for searchableSurface

            //- Range of global indices that can be returned.
            virtual label globalSize() const
            {
                return globalTris().totalSize();
            }

            //- Flip triangles, outsideVolumeType and all cached inside/outside.
            virtual void flip();

            // Other interface bits:
            //
            // regions()
            // hasVolumeType()
            // outsideVolumeType()
            // size()
            // nPoints()
            // coordinates()
            // points()
            // overlaps(const boundBox& bb)

            //- Get bounding spheres (centre and radius squared),
            //- one per element.
            virtual void boundingSpheres
            (
                pointField& centres,
                scalarField& radiusSqr
            ) const;

            virtual void findNearest
            (
                const pointField& sample,
                const scalarField& nearestDistSqr,
                List<pointIndexHit>&
            ) const;

            //- Find the nearest locations for the supplied points to a
            //  particular region in the searchable surface.
            virtual void findNearest
            (
                const pointField& samples,
                const scalarField& nearestDistSqr,
                const labelList& regionIndices,
                List<pointIndexHit>& info
            ) const;

            virtual void findLine
            (
                const pointField& start,
                const pointField& end,
                List<pointIndexHit>&
            ) const;

            virtual void findLineAny
            (
                const pointField& start,
                const pointField& end,
                List<pointIndexHit>&
            ) const;

            //- Get all intersections in order from start to end.
            virtual void findLineAll
            (
                const pointField& start,
                const pointField& end,
                List<List<pointIndexHit>>&
            ) const;

            //- From a set of points and indices get the region
            virtual void getRegion
            (
                const UList<pointIndexHit>&,
                labelList& region
            ) const;

            //- From a set of points and indices get the normal
            virtual void getNormal
            (
                const UList<pointIndexHit>&,
                vectorField& normal
            ) const;

            //- Determine type (inside/outside/mixed) for point. unknown if
            //  cannot be determined (e.g. non-manifold surface)
            virtual void getVolumeType
            (
                const pointField& points,
                List<volumeType>& volType
            ) const;

            //- Set bounds of surface. Bounds currently set as list of
            //  bounding boxes. Will do redistribution of surface to locally
            //  have all triangles overlapping bounds.
            //  Larger bounds: more triangles (memory), more fully local tests
            //  (quick).
            //  keepNonLocal = true : keep triangles that do not overlap
            //  any processor bounds.
            //  Should really be split into a routine to determine decomposition
            //  and one that does actual distribution but determining
            //  decomposition with duplicate triangle merging requires
            //  same amount as work as actual distribution.
            virtual void distribute
            (
                const UList<treeBoundBox>& followBb,
                const bool keepNonLocal,
                autoPtr<mapDistribute>& faceMap,
                autoPtr<mapDistribute>& pointMap
            );


        // Other

            //- WIP. From a set of hits (points and indices)
            //- get the specified field. Misses do not get set.
            virtual void getField
            (
                const UList<pointIndexHit>&,
                labelList&
            ) const;

            //- WIP.
            //- Retrieve the "values" field entries (if it exists)
            virtual label getField(labelList& values) const;

            //- WIP.
            //- Query for surface information.
            //- \returns the number of surface elements queried.
            //- \returns 0 for an empty surface or if all query parameters
            //- were nullptr.
            virtual label measure
            (
                //- [out] surface element centres
                List<vector>* centres,
                //- [out] surface element bounding sphere radius(squared)
                List<scalar>* boundingRadiusSqr,
                //- [out] surface element region ids
                List<label>* regions,
                //- [out] surface element unit normals
                List<vector>* normals
            ) const;

            //- Obtains global indices from pointIndexHit and swaps them back
            //  to their original processor. Used to calculate local region
            //  and normal.
            virtual autoPtr<mapDistribute> localQueries
            (
                const UList<pointIndexHit>& info,
                labelList& triangleIndex
            ) const;

            //- Print some stats.
            //- Parallel aware version of triSurface::writeStats.
            void writeStats(Ostream& os) const;

            //- Print some information about triangles and bounds.
            void printBounds() const;


        // regIOobject implementation

            //- Write using stream options
            //  Do not use the triSurfaceMesh::writeObject since it
            //  would filter out empty regions. These need to be preserved
            //  in case we want to make decisions based on the number of
            //  regions.
            virtual bool writeObject
            (
                IOstreamOption streamOpt,
                const bool writeOnProc
            ) const;

            //- Is object global
            virtual bool global() const
            {
                return false;
            }

            //- Return complete path + object name if the file exists
            //  either in the case/processor or case otherwise null
            virtual fileName filePath() const
            {
                return searchableSurface::localFilePath(type());
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Global file type for distributedTriSurfaceMesh
template<>
struct is_globalIOobject<distributedTriSurfaceMesh> : std::true_type {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "distributedTriSurfaceMesh.txx"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
