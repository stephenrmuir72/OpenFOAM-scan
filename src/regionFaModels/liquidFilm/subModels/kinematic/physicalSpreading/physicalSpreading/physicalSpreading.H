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
    Foam::regionModels::areaSurfaceFilmModels::physicalSpreading

Description
    Abstract base class for conservative physical spreading of
    Lagrangian-to-film mass, momentum and normal-impulse sources.

Usage
    // Minimal example
    \verbatim
    physicalSpreading
    {
        // Mandatory entries
        type                <word>;

        // Optional entries
        beta                <Function1<scalar>>;
        footprintShape      <word>;
        footprintDistance   <word>;
        footprintKernel     <word>;
        faceCoverage        <word>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property        | Description            | Type | Reqd | Deflt
      footprintShape  | Footprint-shape model  | word | no   | circular
      footprintDistance | Footprint-distance mdl | word | no | localTangentPlane
      footprintKernel | Footprint-kernel model | word | no | uniform
      faceCoverage    | Face-coverage model    | word | no | faceCentre
    \endtable

    The inherited entries are elaborated in:
      - \link footprintShape.H \endlink
      - \link footprintDistance.H \endlink
      - \link footprintKernel.H \endlink
      - \link faceCoverage.H \endlink

Note
  - The model is optional. Absence of the dictionary disables spreading.

SourceFiles
    physicalSpreading.C
    physicalSpreadingNew.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_physicalSpreading_H
#define regionFaModels_physicalSpreading_H

#include "areaFields.H"
#include "autoPtr.H"
#include "dictionary.H"
#include "faMesh.H"
#include "physicalSpreadingGeometry.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// Forward Declarations
class footprintShape;
class footprintDistance;
class footprintKernel;
class faceCoverage;

/*---------------------------------------------------------------------------*\
                      Class physicalSpreading Declaration
\*---------------------------------------------------------------------------*/

class physicalSpreading
{
protected:

    // Protected Data

        //- Const reference to the finite-area mesh
        const faMesh& mesh_;

        //- Model dictionary for run-time input diagnostics
        const dictionary dict_;

        //- Shared read-only global finite-area geometry and adjacency,
        //- consumed by all four component models
        physicalSpreadingGeometry geometry_;

        //- Selected footprint-shape component model
        autoPtr<footprintShape> shapePtr_;

        //- Selected footprint-distance component model
        autoPtr<footprintDistance> distancePtr_;

        //- Selected distribution-kernel component model
        autoPtr<footprintKernel> kernelPtr_;

        //- Selected face-coverage component model
        autoPtr<faceCoverage> coveragePtr_;


    // Protected Member Functions

        //- Spread all source channels using one shared set of weights,
        //- delegating shape/distance/coverage/kernel evaluation to the
        //- selected component models
        void spreadSources
        (
            areaScalarField& mass,
            areaVectorField& momentum,
            areaScalarField& normalImpulse,
            const areaScalarField& diameter,
            const scalar beta
        ) const;


public:

    //- Runtime type information
    TypeName("physicalSpreading");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            physicalSpreading,
            dictionary,
            (
                const faMesh& mesh,
                const dictionary& dict
            ),
            (mesh, dict)
        );


    // Constructors

        //- Construct from model type, mesh and dictionary
        physicalSpreading
        (
            const word& modelType,
            const faMesh& mesh,
            const dictionary& dict
        );


    // Selectors

        //- Select a physical-spreading model
        static autoPtr<physicalSpreading> New
        (
            const faMesh& mesh,
            const dictionary& dict
        );

        //- Select the optional physicalSpreading sub-dictionary
        static autoPtr<physicalSpreading> NewIfPresent
        (
            const faMesh& mesh,
            const dictionary& filmDict
        );


    // Generated Methods

        //- No copy construct
        physicalSpreading(const physicalSpreading&) = delete;

        //- No copy assignment
        void operator=(const physicalSpreading&) = delete;


    //- Destructor
    virtual ~physicalSpreading();


    // Member Functions

        //- Spread mass, momentum and normal impulse together
        virtual void spread
        (
            areaScalarField& mass,
            areaVectorField& momentum,
            areaScalarField& normalImpulse,
            const areaScalarField& diameter
        ) const = 0;

        //- Flag indicating whether parcel diameter mapping is required
        virtual bool needsDiameter() const noexcept { return false; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
