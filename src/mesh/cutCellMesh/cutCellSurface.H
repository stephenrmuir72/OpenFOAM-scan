/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
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

Class
    Foam::cutCellSurface

Description
    A triangulated surface to cut a mesh with, built from one entry of the
    geometry dictionary.

    Two kinds of entry are supported in this slice. A searchableBox becomes
    exact triangles, so the S1 test cases stay dictionary-only, and a surface
    file is read as a triSurface. Curved dictionary shapes are generated offline
    as STL (07 section 6.11), so this is the complete set.

    The cutting code reasons about which side of the surface a piece lies on,
    which only means anything for a closed, consistently oriented surface. Both
    are checked at construction and reported rather than left to produce a
    plausible but wrong mesh:

    - every directed edge appears exactly once, and its reverse exactly once,
      which holds for a closed surface and fails for an open one, and
    - the signed volume is positive, which distinguishes outward orientation
      from inward.

SourceFiles
    cutCellSurface.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_cutCellSurface_H
#define Foam_cutCellSurface_H

#include "word.H"
#include "pointField.H"
#include "triFace.H"
#include "triSurface.H"
#include "triSurfaceSearch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dictionary;
class Time;

/*---------------------------------------------------------------------------*\
                      Class cutCellSurface Declaration
\*---------------------------------------------------------------------------*/

class cutCellSurface
{
    // Private Data

        //- Name of the entry in the geometry dictionary
        word name_;

        //- Geometry type from the dictionary
        word type_;

        //- Triangulation, outward oriented
        triSurface surface_;

        //- Search structure over surface_, which refers to it
        triSurfaceSearch search_;

        //- Names of the surface's regions, in the order of their first
        //- appearance in the file. Empty when the format carries no names,
        //- which is the case for an STL and for a dictionary shape.
        wordList regionNames_;

        //- Region of each triangle, indexed as surface_ is. Kept here because
        //- the region is lost twice on the way to surface_: once when the
        //- triangles are read into triFaces, and again by the triSurface
        //- constructor from triFaces, which puts every face in region 0.
        labelList regionOf_;


    // Private Member Functions

        //- Path of the entry's surface file
        static fileName surfacePath
        (
            const word& name,
            const dictionary& dict,
            const Time& runTime
        );

        //- Region names and per-face regions from the file. Both empty if the
        //- format carries neither.
        static void readRegions
        (
            const word& name,
            const dictionary& dict,
            const Time& runTime,
            wordList& names,
            labelList& regionOf
        );

        //- Triangles of a searchableBox entry, outward oriented
        static void boxTriangles
        (
            const word& name,
            const dictionary& dict,
            pointField& points,
            List<triFace>& tris
        );

        //- Triangles read from a surface file
        static void fileTriangles
        (
            const word& name,
            const dictionary& dict,
            const Time& runTime,
            pointField& points,
            List<triFace>& tris
        );

        //- Fail unless the surface is closed and outward oriented
        static void checkClosedOriented(const triSurface& s, const word& name);

        //- Signed volume of a triangulation, referenced to its bounding box
        //- centre so that a surface far from the origin is safe
        static scalar signedVolume(const triSurface& s);

        //- Build the triangulation of a geometry dictionary entry. Checks the
        //- result and reports, so the constructor only has to hold it.
        static triSurface makeSurface
        (
            const word& name,
            const dictionary& dict,
            const Time& runTime
        );


public:

    // Constructors

        //- Construct from a geometry dictionary entry
        cutCellSurface
        (
            const word& name,
            const dictionary& dict,
            const Time& runTime
        );


    // Member Functions

        //- Name of the entry
        const word& name() const noexcept
        {
            return name_;
        }

        //- True when the geometry is a searchableBox
        bool isSearchableBox() const noexcept
        {
            return type_ == "searchableBox";
        }

        //- Triangulation
        const triSurface& surface() const noexcept
        {
            return surface_;
        }

        //- Region names, empty when the surface's format carries none
        const wordList& regionNames() const noexcept
        {
            return regionNames_;
        }

        //- Region of each triangle, indexed as surface() is
        const labelList& regionOf() const noexcept
        {
            return regionOf_;
        }

        //- Search structure
        const triSurfaceSearch& search() const noexcept
        {
            return search_;
        }

        //- Number of triangles
        label size() const noexcept
        {
            return surface_.size();
        }

        //- Points of one triangle
        pointField triangle(const label i) const;

        //- Signed volume, positive when the surface is outward oriented
        scalar signedVolume() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
