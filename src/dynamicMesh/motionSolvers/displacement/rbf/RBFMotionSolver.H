/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2009 Frank Bos, TU Delft
    Copyright (C) 2016-2020 Applied CCM Pty Ltd
    Copyright (C) 2026 Keysight Technologies
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
    Foam::RBFMotionSolver

Description
    Radial basis function motion solver.

Usage
    Example of the dictionary specification:
    \verbatim
    {
        motionSolver    RBFMotionSolver;
        interpolation
        {
            RBF         WendlandC2;
            radius      <scalar>;
            ...
        }
        movingPatches <patch selection>;
        staticPatches <patch selection>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property      | Description                       | Type   | Reqd | Deflt
      motionSolver  | type: RBFMotionSolver             | word   | yes  | -
      interpolation | RBF control parameters            | dict   | yes  | -
      movingPatches | Selection of motion patches       | wordRes | yes | -
      staticPatches | Selection of zero-motion patches  | wordRes | yes | -
      coarseningRatio | Control point coarsening        | label  | yes  | -
      frozenInterpolation | No update interpolation parameters | bool | yes  | -
      includeStaticPatches | Include zero-motion in interpolation | bool | yes  | -
      functionObject | Motion driven by a function object | bool | no | false
    \endtable

SeeAlso
    Foam::RBFFunction
    Foam::RBFInterpolation

SourceFiles
    RBFMotionSolver.cxx

\*---------------------------------------------------------------------------*/

#ifndef Foam_RBFMotionSolver_H
#define Foam_RBFMotionSolver_H

#include "points0MotionSolver.H"
#include "IndirectList.H"
#include "RBFInterpolation.H"
#include "solidBodyMotionFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class globalIndex;

/*---------------------------------------------------------------------------*\
                      Class RBFMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class RBFMotionSolver
:
    public points0MotionSolver
{
    // Private Data

        //- Motion function
        autoPtr<solidBodyMotionFunction> SBMFPtr_;

        //- The moving patches selector
        wordRes movingPatches_;

        //- The static patches selector
        wordRes staticPatches_;

        //- Coarsening ratio
        label coarseningRatio_;

        //- Include static patches (zero displacement) in RBF interpolation
        bool includeStaticPatches_;

        //- Frozen interpolation
        bool frozenInterpolation_;

        //- Motion is updated from functionObject (default: false)
        bool useFunctionObject_;

        //- The mesh point IDs from the moving patches [rank local]
        labelList movingIDs_;

        //- The mesh point IDs from the static patches (zero motion)
        //- [rank local]
        labelList staticIDs_;

        //- The point IDs of the interpolated points [rank local].
        //- All zoneMotion pointIDs, excluding moving and static points.
        labelList interpIDs_;

        //- Control point IDs [rank local].
        //- A coarsened subset of movingIDs_ and (optionally) staticIDs_.
        labelList controlIDs_;

        //- The rank local unique contribution to globalControlPoints_
        //- after topological merge
        //- (parallel only).
        labelList uniqControlIDs_;

        //- The global numbering for the coordination between uniqControlIDs_
        //- and globalControlPoints_ (parallel only)
        autoPtr<globalIndex> ctrlPointsNumberingPtr_;

        //- All boundary control points [multiple ranks]
        pointField globalControlPoints_;

        //- The points for interpolation [rank local]
        pointField interpPoints_;

        //- Prescribed motion for points on the moving patches,
        //- corresponding to movingIDs_ locations [rank local]
        vectorField motion_;

        //- The RBF interpolation: holds references to globalControlPoints_
        //- and interpPoints_
        RBFInterpolation interpolation_;


    // Private Member Functions

        //- Set the control points (serial or parallel)
        void setControlPoints();

        //- Make control point IDs.  Constructor helper
        void makeControlIDs();


public:

    //- Runtime type information
    TypeName("RBFMotionSolver");


    // Generated Methods

        //- No copy construct
        RBFMotionSolver(const RBFMotionSolver&) = delete;

        //- No copy assignment
        void operator=(const RBFMotionSolver&) = delete;


    // Constructors

        //- Construct from polyMesh
        RBFMotionSolver(const polyMesh& mesh, const IOdictionary& dict);


    //- Destructor
    virtual ~RBFMotionSolver();

    // Member Functions

        //- The indices of the polyMesh points on moving boundaries
        const labelUList& movingIDs() const noexcept { return movingIDs_; }

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Update local data for geometry changes
        virtual void movePoints(const pointField&);

        //- Update the mesh corresponding to given map
        virtual void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
