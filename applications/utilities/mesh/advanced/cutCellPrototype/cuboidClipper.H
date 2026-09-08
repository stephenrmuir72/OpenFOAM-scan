/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O p e r a t i o n |
    \\  /    A n d           |
     \\/     M a n i p u l a t i o n  |
-------------------------------------------------------------------------------
    Description
        Utility-owned clipping of a cut-cell interface for an axis-aligned
        six-region cuboid.
\*---------------------------------------------------------------------------*/

#ifndef cuboidClipper_H
#define cuboidClipper_H

#include "Map.H"
#include "faceList.H"
#include "labelList.H"
#include "scalar.H"

namespace Foam
{

class meshCutter;
class mapPolyMesh;
class polyMesh;
class triSurface;

/*---------------------------------------------------------------------------*\
                       Class cuboidClipper Declaration
\*---------------------------------------------------------------------------*/

class cuboidClipper
{
public:

    //- Rebuild an axis-aligned support mesh as coherent cuboid intersections.
    static void applyDirect
    (
        polyMesh& mesh,
        const triSurface& surface,
        const scalar tolerance,
        const label keepZone
    );

    //- Replace multi-region meshCutter faces with planar cuboid-region faces.
    //  faceRegions is updated to the post-change face labels and contains
    //  provenance for all generated cut-interface faces.
    static void apply
    (
        polyMesh& mesh,
        const triSurface& surface,
        const meshCutter& cutter,
        const labelList& originalCutCells,
        const labelList& mappedCutCells,
        const List<labelList>& originalCellFaces,
        const List<faceList>& originalCellFaceGeometry,
        const mapPolyMesh& cutMap,
        const List<labelList>& cutCellRegions,
        const scalar tolerance,
        const bool retainedOnly,
        Map<label>& faceRegions,
        Map<label>& retainedCells
    );
};

} // End namespace Foam

#endif

// ************************************************************************* //
