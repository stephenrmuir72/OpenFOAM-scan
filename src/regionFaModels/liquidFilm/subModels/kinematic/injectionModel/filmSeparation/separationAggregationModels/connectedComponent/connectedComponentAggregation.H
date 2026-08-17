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
    Foam::regionModels::areaSurfaceFilmModels::connectedComponentAggregation

Description
    Conservative aggregation of finite-area film-separation sources over
    connected active face components.

    The model decouples the film stripping rate from the Lagrangian parcel
    injection frequency by buffering separated mass and momentum per face.
    A Lagrangian parcel is injected when a physically coherent connected
    component reaches its mass-weighted event time. This prevents the
    "particle explosion" problem where
    \f$N_{\text{parcels}} \propto N_{\text{faces}} \times N_{\text{steps}}\f$.

    The model collects buffered separated mass and momentum over physically
    connected finite-area faces, selects a deterministic representative face
    closest to the mass centroid, injects the aggregate source there, and
    resets all member-face accumulators.  Processor-boundary connections are
    included via finite-area boundary connectivity and global face numbering.
    Two neighbouring faces are joined only when their morphology agrees and
    their centre separation is within both local physical coherence reaches,
    so the grouping is independent of decomposition and dimensional mesh
    thresholds.

Usage
    Minimal example:
    \verbatim
    separationAggregation
    {
        type                   connectedComponent;

        scaleModel
        {
            type morphologyBased;
            ...
        }

        diameterModel          <word>;
        fixedDiameter          <scalar>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property      | Description                      | Type   | Reqd | Deflt
      type          | Type name: connectedComponent    | word   | yes  | -
      diameterModel | Aggregate parcel diameter model               <!--
                    -->             | word | no | massWeightedFilmThickness
      fixedDiameter | Fixed parcel diameter [m]        | scalar | no | 0
    \endtable

SourceFiles
    connectedComponentAggregation.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_connectedComponentAggregation_H
#define regionFaModels_connectedComponentAggregation_H

#include "separationAggregationModel.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                 Class connectedComponentAggregation Declaration
\*---------------------------------------------------------------------------*/

class connectedComponentAggregation
:
    public separationAggregationModel
{
public:

    // Public Enumerations

        //- Options for the aggregate parcel diameter
        enum diameterModelType : char
        {
            MASS_WEIGHTED_FILM_THICKNESS = 0,  //!< Use sum(m*d)/sum(m)
            FIXED_DIAMETER                     //!< Use a prescribed diameter
        };

        //- Names for diameterModelType
        static const Enum<diameterModelType> diameterModelTypeNames;


private:

    // Private Data

        //- Diameter model
        diameterModelType diameterModel_;

        //- Fixed parcel diameter [m]
        scalar fixedDiameter_;


public:

    //- Runtime type information
    TypeName("connectedComponent");


    // Constructors

        //- Construct from film and dictionary
        connectedComponentAggregation
        (
            const liquidFilmBase& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~connectedComponentAggregation() = default;


    // Member Functions

        //- True for the legacy fixed aggregate-diameter option
        virtual bool overridesDiameter() const noexcept
        {
            return diameterModel_ == FIXED_DIAMETER;
        }

        //- Aggregate buffered separated mass into injection fields
        virtual void aggregate
        (
            scalarField& massAccumulator,
            vectorField& momentumAccumulator,
            scalarField& diameterAccumulator,
            scalarField& accumulationTime,
            const scalarField& eventTime,
            const scalarField& coherenceLength,
            const labelField& morphology,
            scalarField& massToInject,
            scalarField& diameterToInject,
            vectorField& velocityToInject
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
