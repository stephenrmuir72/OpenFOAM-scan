/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::cornerDetectionModels::geometryBased

Description
    Geometry-based corner detection model. Identifies sharp or round edges from
    face-normal dihedral angles and curvature, then marks separation faces.

    With the default \c fluxDirected face selection, only the upstream face at
    each corner edge is flagged using the film edge flux \c phi2s (same rule as
    \c fluxBased). The legacy \c symmetric selection marks both faces sharing a
    sharp edge.

Usage
    Minimal example:
    \verbatim
    filmSeparationCoeffs
    {
        // Inherited entries
        ...

        // Optional entries
        cornerDetectionModel  geometryBased;
        cornerCurveType       <word>;
        cornerType            <word>;
        angleSharp            <scalar>;   // [deg]
        angleRoundMin         <scalar>;   // [deg]
        angleRoundMax         <scalar>;   // [deg]
        maxRoundRadius        <scalar>;   // [m]
        minEdgeLength         <scalar>;   // [m]
        nSmooth               <label>;    // [no. of passes]
        sharpBoundaryEdges    <bool>;
        faceSelection         <word>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                   | Type | Reqd | Deflt
      cornerDetectionModel | Corner detector model | word | no   | fluxBased
      cornerCurveType | Corner-curvature type      | word | no   | any
      cornerType      | Corner type                | word | no   | sharpOrRound
      angleSharp      | Sharp-angle limit [deg]    | scalar | no | 45
      angleRoundMin   | Minimum round-angle limit [deg] | scalar | no | 5
      angleRoundMax   | Maximum round-angle limit [deg] | scalar | no | 45
      maxRoundRadius  | Maximum round-radius limit [m] | scalar | no | 2e-3
      minEdgeLength   | Minimum edge length [m]     | scalar | no | 0
      nSmooth         | No. of smoothing passes     | label | no | 0
      sharpBoundaryEdges | Treat boundary edges as sharp | bool | no | false
      faceSelection   | Corner-face marking strategy | word | no | fluxDirected
    \endtable

    Options for the \c faceSelection entry:
    \verbatim
        fluxDirected | Mark only the upstream separation face (phi2s sign)
        symmetric    | Mark both faces sharing a sharp edge (legacy/debug)
    \endverbatim

    Options for the \c cornerCurve entry:
    \verbatim
        any          | Convex or concave corners
        convex       | Convex corners
        concave      | Concave corners
    \endverbatim

    Options for the \c cornerType entry:
    \verbatim
        sharp        | Sharp corners
        round        | Round corners
        sharpOrRound | Sharp or round corners
    \endverbatim

    The inherited entries are elaborated in:
      - \link cornerDetectionModel.H \endlink

SourceFiles
    geometryBased.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_cornerDetectionModels_geometryBased_H
#define Foam_cornerDetectionModels_geometryBased_H

#include "cornerDetectionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace cornerDetectionModels
{

/*---------------------------------------------------------------------------*\
                          Class geometryBased Declaration
\*---------------------------------------------------------------------------*/

class geometryBased
:
    public cornerDetectionModel
{
    // Private Enumerations

        //- Options for the corner-curvature type
        enum cornerCurveType : char
        {
            ANY = 0,    //!< "Convex or concave corners"
            CONVEX,     //!< "Convex corners"
            CONCAVE     //!< "Concave corners"
        };

        //- Names for cornerCurveType
        static const Enum<cornerCurveType> cornerCurveTypeNames;


        //- Options for the corner type
        enum cornerType : char
        {
            ALL = 0,    //!< "Sharp or round corners"
            SHARP,      //!< "Sharp corners"
            ROUND       //!< "Round corners"
        };

        //- Names for cornerType
        static const Enum<cornerType> cornerTypeNames;


        //- Options for corner-face marking
        enum faceSelectionType : char
        {
            FLUX_DIRECTED = 0,  //!< "fluxDirected"
            SYMMETRIC              //!< "symmetric"
        };

        //- Names for faceSelectionType
        static const Enum<faceSelectionType> faceSelectionTypeNames_;


    // Private Data

        //- Corner-curvature type
        enum cornerCurveType cornerCurveType_;

        //- Corner type
        enum cornerType cornerType_;

        //- Flag to deduce if the object is initialised
        bool init_;

        //- Bitset of edges identified as sharp
        bitSet sharpEdges_;

        //- Bitset of edges identified as round
        bitSet roundEdges_;

        //- Bitset of edges identified as a combination of sharp and round
        bitSet cornerEdges_;

        //- Bitset of faces identified as sharp
        bitSet sharpFaces_;

        //- Bitset of faces identified as round
        bitSet roundFaces_;

        //- Bitset of faces identified as a combination of sharp and round
        bitSet cornerFaces_;

        //- Sharp-angle limit
        scalar angleSharpDeg_;

        //- Minimum round-angle limit
        scalar angleRoundMinDeg_;

        //- Maximum round-angle limit
        scalar angleRoundMaxDeg_;

        //- Maximum round-radius limit
        scalar maxRoundRadius_;

        //- Minimum edge length; ignore edges shorter than this
        scalar minEdgeLength_;

        //- Number of smoothing passes on the binary edge mask
        label nSmooth_;

        //- Flag to treat one-face boundary edges as sharp
        bool sharpBoundaryEdges_;

        //- Corner-face marking strategy
        enum faceSelectionType faceSelection_;


    // Private Member Functions

        //- Return the signed convex(+1)/concave(-1) bending sense across an
        //- edge using the decomposition-invariant face-centre form
        //- sign((fcN - fcO) & (n1 - n0)) (symmetric under owner/neighbour swap)
        scalar curvatureSign
        (
            const vector& fcO,
            const vector& fcN,
            const vector& n0,
            const vector& n1
        ) const;

        //- Classify edges into sharp/round according to dihedral angle and
        //- inferred radius
        void classifyEdges
        (
            bitSet& sharpEdges,
            bitSet& roundEdges
        ) const;

        //- Convert an edge mask to a face mask (face is set if any of its
        //- edges are set)
        void edgesToFaces
        (
            const bitSet& edgeMask,
            bitSet& faceMask
        ) const;

        //- Mark separation faces from corner edges using film edge flux
        bitSet edgesToSeparationFaces(const bitSet& edgeMask) const;

        //- Return the list of corner angles [rad] for each edge
        scalarList calcCornerAngles
        (
            const bitSet& faces,
            const bitSet& edges
        ) const;

        // Write edges and faces as OBJ sets for debug purposes
        void writeEdgesAndFaces() const;


public:

    //- Runtime type information
    TypeName("geometryBased");


    // Constructors

        //- Construct from components
        geometryBased
        (
            const faMesh& mesh,
            const regionModels::areaSurfaceFilmModels::liquidFilmBase& film,
            const dictionary& dict
        );


    // Destructor
    virtual ~geometryBased() = default;


    // Member Functions

    // Evaluation

        //- Detect and store the corner faces
        virtual bool detectCorners();


    // I-O

        //- Read the model dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace cornerDetectionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
