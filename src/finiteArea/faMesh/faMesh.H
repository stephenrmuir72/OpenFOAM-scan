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

Class
    Foam::faMesh

Description
    Finite area mesh (used for 2-D non-Euclidian finite area method).

    In the finiteArea refactor an faMesh is a reduced-order (2D) manifold
    fvMesh whose "cells" are the surface polygons and whose "faces" are the
    bounding edges. Being a thin fvMesh subclass it reuses the entire
    finite-volume container/matrix/interpolation/boundary/parallel machinery
    directly; the finite-area schemes operate on it through the fvMesh
    interface.

    In addition it owns the (optional) coupling addressing that links the
    area mesh back to the boundary faces of a parent volume mesh. This
    provides a stable home for the mapping used by volSurfaceMapping and
    film/shell applications, which the pure fvMesh does not have.

SourceFiles
    faMesh.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_faMesh_H
#define Foam_faMesh_H

#include "fvMesh.H"
#include "labelList.H"
#include "labelPair.H"
#include "boolList.H"
#include "vectorField.H"
#include "scalarField.H"
#include "uindirectPrimitivePatch.H"
#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyBoundaryMesh;
class wordRes;

/*---------------------------------------------------------------------------*\
                           Class faMesh Declaration
\*---------------------------------------------------------------------------*/

class faMesh
:
    public fvMesh
{
    // Private Data

        //- The area-region name
        word areaName_;

        //- Parent-mesh face labels that this area mesh maps onto (optional;
        //  empty for a standalone area region)
        labelList meshFaces_;

        //- Parent-mesh point label for each area-mesh point (optional; empty
        //  for a standalone area region). Used by the parallel addressing shim
        //  to derive pointProcAddressing from the parent decomposition.
        labelList meshPointLabels_;


    // Demand-driven coupling data

        //- The polyPatch/local-face for each of the meshFaces()
        mutable std::unique_ptr<List<labelPair>> polyPatchFacesPtr_;

        //- Sorted list of parent polyPatch ids referenced by meshFaces()
        mutable std::unique_ptr<labelList> polyPatchIdsPtr_;

        //- Demand-driven area patch (reduced-order face polygons) on the
        //- parent mesh, addressed by meshFaces()
        mutable std::unique_ptr<uindirectPrimitivePatch> patchPtr_;

        //- Per finite-area boundary edge (reduced-order boundary face), the
        //  parent-mesh neighbour boundary face across that edge (-1 if none,
        //  e.g. an inter-processor edge whose neighbour area face is on another
        //  rank). Indexed by reduced-order boundary-face index
        //  (facei - nInternalFaces()). Drives the ngbPolyPatch* accessors.
        mutable std::unique_ptr<labelList> ngbPolyPatchFacesPtr_;

        //- Point-normal-correction flag per boundary patch (finiteArea compat
        //  hook, see correctPatchPointNormals())
        mutable std::unique_ptr<boolList> correctPatchPointNormalsPtr_;

        //- Per parent finite-area patch edge (patch().edges() index), the
        //  reduced-order boundary face (global mesh-face id) coincident with
        //  it, or -1 for an internal edge. Bridges the parent-patch edge
        //  topology to the reduced-order boundary patches.
        mutable std::unique_ptr<labelList> patchEdgeToBndFacePtr_;

        //- Per reduced-order point, the corresponding patch() local point
        //  index. The reduced-order mesh and patch() enumerate the same
        //  free-surface points (both keyed on parent-mesh points) but not
        //  necessarily in the same order, so this map is required whenever a
        //  reduced-order face's point labels index a patch-local field
        //  (points, oldPoints, ...). Identity for a standalone/consistent
        //  region.
        mutable std::unique_ptr<labelList> patchPointMapPtr_;


    // Private Static Data

        //- The prefix to the parent registry name: "finite-area"
        static const word prefix_;


    // Private Member Functions

        //- Calculate the whichPatchFaces()/whichPolyPatches() addressing
        void calcWhichPatchFaces() const;

        //- Calculate ngbPolyPatchFacesPtr_: for each finite-area boundary edge
        //  the parent-mesh neighbour boundary face across it. The neighbour of
        //  a genuine (e.g. wall) rim edge is a boundary face of the SAME parent
        //  domain (the area faces are boundary faces of the local domain, and a
        //  wall face sharing a rim edge is local too), so no off-processor
        //  exchange is needed. No-op (all -1) for a standalone area region.
        void calcNgbPolyPatchFaces() const;

        //- Calculate patchEdgeToBndFacePtr_ (parent-patch edge -> reduced-order
        //  boundary face)
        void calcPatchEdgeToBoundaryFace() const;

        //- Calculate patchPointMapPtr_ (reduced-order point -> patch() local
        //  point index)
        void calcPatchPointMap() const;

        //- Publish, in this mesh's registry, the parent-mesh point label of
        //  each reduced-order point (a labelIOList named
        //  "pointAreaNormalsParentMap"), so the reduced-order (2D) geometry
        //  scheme can sum point-area-normal contributions through the PARENT
        //  mesh's correct global point coupling instead of the reduced-order
        //  mesh's own (degenerate) coupling. See
        //  basicFvGeometryScheme::calcPointAreaNormals for why the reduced-order
        //  mesh's own coupling corrupts the normal at multiply-shared points
        //  (e.g. the poles of a decomposed sphere). No-op when standalone (no
        //  parent) or when the parent point addressing is unavailable.
        void registerParentPointMap();

        //- (Re)publish, in this mesh's registry, the per-reduced-order-point
        //  neighbour (wall) polyPatch face normal to be removed from the
        //  point-area normal at contact-line points (a vectorIOField named
        //  "pointAreaNormalsBoundaryCorrection"), for patches opted in via
        //  correctPatchPointNormals(). The reduced-order pointAreaNormals scheme
        //  (basicFvGeometryScheme::calcPointAreaNormals, finiteVolume) applies it
        //  inside the point-normal build so it reaches faceCurvatures, mirroring
        //  the legacy faMesh::calcPointAreaNormals haloFaceNormals correction.
        //  No-op when standalone or when no patch requests the correction.
        void updatePointNormalsBoundaryCorrection() const;

        //- Finalise the reduced-order boundary in a parallel-safe order.
        //  Drives the boundary topology/geometry collectives explicitly and
        //  uniformly on every rank, so an empty area subdomain issues the
        //  same worldComm collectives as the ranks that carry the region
        //  (see the implementation for the deadlock this avoids).
        void finaliseParallelBoundary();

        //- Repopulate meshFaces_/meshPointLabels_ from the persisted
        //- finite-area[/<region>]/faMesh/{faceLabels,pointLabels} written by
        //- writeObject(). Needed by the IOobject (read) constructor - whose
        //  fvMesh base only reads the reduced-order polyMesh, not the parent
        //  coupling addressing - so that faceLabels()/pointLabels()/patch()
        //  (and hence volSurfaceMapping / faMeshProcAddressing) work after a
        //  createFaMesh-style read. No-op when standalone (no parent) or when
        //  the addressing is already populated. pointLabels is read
        //  READ_IF_PRESENT so meshes written before it was persisted still
        //  load (with an empty meshPointLabels_, as before).
        void readCouplingAddressing();

        //- Bundle used to build the fvMesh base + coupling in one ctor
        struct meshData
        {
            autoPtr<fvMesh> region;
            labelList faceLabels;

            //- Parent-mesh point label for each area-mesh point
            labelList pointLabels;

            //- Optional definition-derived boundary patch names/types
            //  (ordered by the region boundary patch index). When non-empty
            //  the delegated ctor builds named/typed patches instead of
            //  cloning the parent-derived region patches.
            wordList patchNames;
            wordList patchTypes;

            //- Neighbour processor for each emitted patch (-1 for a
            //  non-processor patch). For patchTypes[i] == "processor" this is
            //  the rank on the other side of the inter-domain edges, used to
            //  build a genuine fvMesh processorPolyPatch on the reduced-order
            //  region.
            labelList patchNeighbProcs;
        };

        //- Read the parent-mesh faceLabels for the given area region
        static labelList readFaceLabels
        (
            const word& areaName,
            const polyMesh& pMesh
        );

        //- Select the parent boundary faceLabels matching the given patches
        static labelList selectPatchFaces
        (
            const polyBoundaryMesh& pbm,
            const wordRes& patchNames
        );

        //- Build a temporary reduced-order fvMesh from the parent mesh and
        //- the selected parent faceLabels
        static autoPtr<fvMesh> buildRegionMesh
        (
            const word& areaName,
            const polyMesh& pMesh,
            const labelList& faceLabels,
            labelList& pointLabels
        );

        //- Build a temporary reduced-order fvMesh, grouping the boundary
        //- faces (edges) into the named/typed patches declared in the
        //- faMeshDefinition "boundary" sub-dictionary. The resulting patch
        //- names/types (ordered by region boundary patch index) are returned
        //- via \p patchNames / \p patchTypes.
        static autoPtr<fvMesh> buildRegionMeshFromDefinition
        (
            const word& areaName,
            const polyMesh& pMesh,
            const labelList& faceLabels,
            const dictionary& faMeshDefinition,
            wordList& patchNames,
            wordList& patchTypes,
            labelList& patchNeighbProcs,
            labelList& pointLabels
        );

        //- Read faceLabels + build the reduced-order region mesh
        static meshData buildMeshData
        (
            const word& areaName,
            const polyMesh& pMesh
        );

        //- Generate faceLabels + reduced-order region mesh from a definition
        static meshData buildMeshDataFromDefinition
        (
            const word& areaName,
            const polyMesh& pMesh,
            const dictionary& faMeshDefinition
        );

        //- The IOobject used to register the area region on the parent
        static IOobject regionIO
        (
            const word& areaName,
            const polyMesh& pMesh
        );

        //- Read a legacy finiteArea control dictionary
        //- (system/finite-area/<faDictName>) as fallback content for the
        //- fvSchemes/fvSolution controls. Returns an empty dictionary when
        //- the file is absent. Used to preserve backward compatibility with
        //- cases that ship faSchemes/faSolution (D2).
        static dictionary readLegacyControlDict
        (
            const word& faDictName,
            const objectRegistry& db
        );

        //- Delegated constructor: build fvMesh base + record coupling
        faMesh
        (
            const word& areaName,
            const polyMesh& pMesh,
            meshData&& data
        );

        //- No copy construct
        faMesh(const faMesh&) = delete;

        //- No copy assignment
        void operator=(const faMesh&) = delete;


public:

    // Public Typedefs

        //- The mesh type
        typedef faMesh Mesh;


    // Static Data

        //- The mesh sub-directory name (usually "faMesh")
        static word meshSubDir;


    // Constructors

        //- Construct from IOobject (standalone area region, no parent).
        //  Reads the region fvSchemes/fvSolution if present, otherwise falls
        //  back to legacy system/finite-area/faSchemes|faSolution (D2 compat).
        explicit faMesh(const IOobject& io, const bool doInit = true);

        //- Construct from IOobject, capturing the parent-mesh coupling faces
        faMesh
        (
            const IOobject& io,
            const labelUList& meshFaces,
            const bool doInit = true
        )
        :
            fvMesh(io, doInit),
            areaName_(io.name()),
            meshFaces_(meshFaces)
        {}

        //- Construct as a reduced-order region on a parent volume mesh,
        //- reading the parent-mesh faceLabels for the given area region
        faMesh
        (
            const word& areaName,
            const polyMesh& pMesh,
            const bool doInit = true
        );

        //- Construct as a reduced-order region on a parent volume mesh from
        //- a faMeshDefinition dictionary (the makeFaMesh generation path).
        //  Selects the area faces from "polyMeshPatches" and groups the
        //  boundary edges into the named patches under "boundary".
        faMesh
        (
            const word& areaName,
            const polyMesh& pMesh,
            const dictionary& faMeshDefinition,
            const bool doInit = true
        );

        //- Construct the default (single) area region on a parent volume mesh
        //- from a faMeshDefinition dictionary. Convenience overload that uses
        //- the default area-region name (compat with the legacy 2-arg form).
        faMesh
        (
            const polyMesh& pMesh,
            const dictionary& faMeshDefinition,
            const bool doInit = true
        )
        :
            faMesh(word::null, pMesh, faMeshDefinition, doInit)
        {}


    //- Destructor
    virtual ~faMesh() = default;


    // Static Member Functions

        //- The prefix to the parent registry name: "finite-area"
        static const word& prefix() noexcept { return prefix_; }

        //- Return the (optional) parent registry that contains all
        //- finite-area objects on the polyMesh, or nullptr
        static const objectRegistry* registry(const polyMesh& pMesh);

        //- Return the singleton parent registry (on the polyMesh) that
        //- contains all objects related to finite-area
        static const objectRegistry& Registry(const polyMesh& pMesh);

        //- The single-region or specified finite-area region on the polyMesh
        static const faMesh& mesh
        (
            const polyMesh& pMesh,
            const word& areaRegion = word::null
        );

        //- Read-construct an existing area region from the parent polyMesh
        //- when its on-disk faceLabels are available, otherwise return
        //- nullptr. Mirrors the legacy finiteArea faMesh::TryNew selector.
        static autoPtr<faMesh> TryNew
        (
            const word& areaName,
            const polyMesh& pMesh
        );

        //- Read-construct the default (single) area region if available,
        //- otherwise return nullptr.
        static autoPtr<faMesh> TryNew(const polyMesh& pMesh);

        //- True if the on-disk faMesh addressing (faceLabels) for the given
        //- area region exists on the parent polyMesh. Collective: the result
        //- is reduced so all ranks agree.
        static bool hasMeshFiles(const word& areaName, const polyMesh& pMesh);


    // Member Functions

        //- The area-region name
        const word& areaName() const noexcept { return areaName_; }

        //- True if this area mesh is coupled to a parent volume mesh.
        //  Uses the inherited fvMesh::parentMesh() registry lookup.
        bool hasParent() const { return parentMesh().good(); }

        // Keep the inherited fvMesh::movePoints(const pointField&) overload
        // visible alongside the no-argument parent-tracking overload below.
        using fvMesh::movePoints;

        //- Re-slice the reduced-order point positions from the (moved) parent
        //  boundary and propagate the motion to the base fvMesh, clearing and
        //  lazily rebuilding the reduced-order geometry (areaCentres,
        //  faceAreaNormals, pointAreaNormals, faceCurvatures, Le, ...).
        //
        //  The legacy finiteArea faMesh was a parent-polyMesh MoveableMeshObject
        //  and received this automatically from polyMesh::movePoints(). The
        //  finiteArea faMesh is a standalone (reduced-order) fvMesh, so a
        //  moving-mesh consumer (e.g. the free-surface tracker) must call this
        //  after the parent points move; otherwise the reduced-order geometry
        //  and, critically, faceCurvatures() stay frozen at the construction
        //  (flat) configuration. Returns false (no-op) for a standalone region
        //  with no usable parent point map.
        bool movePoints();

        //- The parent-mesh face labels this area mesh maps onto
        const labelList& meshFaces() const noexcept { return meshFaces_; }

        //- The parent-mesh face labels this area mesh maps onto
        //  (compatibility alias for meshFaces())
        const labelList& faceLabels() const noexcept { return meshFaces_; }

        //- The area-mesh face (reduced-order cell) index corresponding to the
        //  given parent-mesh face, or -1 if that parent face is not part of
        //  this area mesh. Inverse lookup into meshFaces()/faceLabels().
        label whichFace(const label meshFacei) const
        {
            return meshFaces_.find(meshFacei);
        }

        //- The parent-mesh point label for each area-mesh point (may be empty
        //  for a standalone area region)
        const labelList& pointLabels() const noexcept
        {
            return meshPointLabels_;
        }

        //- Set the parent-mesh coupling face labels
        void meshFaces(const labelUList& mf) { meshFaces_ = mf; }

        //- The polyPatch/local-face for each of the meshFaces()
        const List<labelPair>& whichPatchFaces() const;

        //- The sorted parent polyPatch ids referenced by meshFaces()
        const labelList& whichPolyPatches() const;

        //- The area patch (reduced-order face polygons) on the parent mesh,
        //  addressed by meshFaces(). Requires a coupled (parent) area mesh.
        const uindirectPrimitivePatch& patch() const;


    // Finite-area boundary <-> parent-polyPatch coupling
    //
    //  The finiteArea reduced-order boundary patches are genuine fvPatches
    //  (processor/empty/wedge coupling is native FV). These accessors restore
    //  the legacy finiteArea faPatch->parent-3D-patch coupling that the
    //  free-surface tracker needs, derived from the parent-mesh addressing this
    //  faMesh already owns rather than a separate boundary halo.

        //- Per finite-area boundary edge (reduced-order boundary face), the
        //  parent-mesh neighbour boundary face across it (-1 if none). Indexed
        //  by reduced-order boundary-face index (facei - nInternalFaces()).
        const labelList& ngbPolyPatchFaces() const;

        //- Per parent finite-area patch edge (patch().edges() index), the
        //  reduced-order boundary face (global mesh-face id) coincident with
        //  it, or -1 for an internal edge. Lets consumers relate the
        //  self-consistent parent-patch edge topology (pointEdges/edgeFaces)
        //  to the reduced-order boundary patches and their coupling data.
        const labelList& patchEdgeToBoundaryFace() const;

        //- Per reduced-order point, the corresponding patch() local point
        //  index. Required to translate a reduced-order face's point labels
        //  (e.g. aMesh().faces()) into the patch-local numbering used by
        //  patch().localPoints()/meshPoints() and fields built on them.
        const labelList& patchPointMap() const;

        //- The parent polyPatch neighbouring the given finite-area boundary
        //  patch (the 3D patch its rim edges sit on), or -1 (e.g. a processor
        //  patch, which has no parent-patch neighbour).
        label ngbPolyPatchIndex(const label faPatchi) const;

        //- Unit normals of the parent faces neighbouring the given finite-area
        //  boundary patch, one per boundary edge (patch face). Zero where there
        //  is no neighbour. Taken directly from the parent faces (local), the
        //  parent-geometry analogue of the legacy faMesh boundary halo.
        tmp<vectorField> ngbPolyPatchFaceNormals(const label faPatchi) const;

        //- Neighbour parent-patch point normals for the given finite-area
        //  boundary patch, one per patch point, summed over the rim edges
        //  incident to each point (mirrors legacy faPatch::ngbPolyPatchPointNormals).
        tmp<vectorField> ngbPolyPatchPointNormals(const label faPatchi) const;

        //- Straight-line point-to-neighbour (P-N) distance per internal
        //  finite-area edge (reduced-order internal face). Used for the
        //  capillary Courant number estimate.
        tmp<scalarField> lPN() const;

        //- The finite-area (reduced-order) points not on the area-mesh boundary
        labelList internalPoints() const;

        //- The finite-area (reduced-order) points on the area-mesh boundary
        labelList boundaryPoints() const;

        //- Point-normal-correction flags per boundary patch (finiteArea compat).
        //  Stored on the mesh so callers (e.g. interfaceTrackingFvMesh's
        //  pointNormalsCorrectionPatches) can request per-patch correction.
        //  Consumed by the reduced-order pointAreaNormals scheme: for a flagged,
        //  non-coupled patch, updatePointNormalsBoundaryCorrection() publishes
        //  the neighbour (wall) polyPatch normal per contact-line point and
        //  basicFvGeometryScheme::calcPointAreaNormals removes it (removeCollinear)
        //  inside the point-normal build, so the finite-area/neighbour-polyPatch
        //  perpendicularity feeds curvature (the legacy haloFaceNormals correction).
        boolList& correctPatchPointNormals() const;

        //- Whether point-normal correction was requested for the given patch
        bool correctPatchPointNormals(const label patchID) const;


    // Finite-area geometry aliases (area == reduced-order cell "volume")

        //- Face areas of the area mesh. On the reduced-order fvMesh the polygon
        //  "cell volume" IS the finite-area face area, so S() aliases V().
        const DimensionedField<scalar, volMesh>& S() const { return V(); }

        //- Old-time face areas (alias for the reduced-order V0())
        const DimensionedField<scalar, volMesh>& S0() const { return V0(); }

        //- Old-old-time face areas (alias for the reduced-order V00())
        const DimensionedField<scalar, volMesh>& S00() const { return V00(); }


    // Finite-area compatibility (edge <-> reduced-order face) addressing

        //- The owner (lower) addressing for the edges (reduced-order faces)
        const labelUList& edgeOwner() const { return owner(); }

        //- The neighbour (upper) addressing for the edges (reduced-order faces)
        const labelUList& edgeNeighbour() const { return neighbour(); }


    // Writing

        //- Write the reduced-order mesh and, for a coupled (parent) area
        //- mesh, the parent-mesh faceLabels addressing
        virtual bool writeObject
        (
            IOstreamOption streamOpt,
            const bool writeOnProc = true
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
