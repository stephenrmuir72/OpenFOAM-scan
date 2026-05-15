/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 OpenCFD Ltd.
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
    Foam::regionModels::areaSurfaceFilmModels::gaussianSpreading

Description
    Spreads Lagrangian-to-film source terms using a Gaussian distribution
    over neighbouring finite-area faces.

    For each face carrying a non-zero source, the stencil is built by
    topological breadth-first search (BFS) for \c nLayers rings of
    neighbours via \c PrimitivePatch::faceFaces() adjacency.  The weight
    for each stencil face \c j relative to the source face \c i is:

    \f[
        w_j = S_j \, \exp\!\left(-\frac{d_{ij}^2}{2\,\sigma^2}\right)
    \f]

    where \f$S_j\f$ is the face area and \f$d_{ij}\f$ the centroid-to-centroid
    distance.  Weights are normalised so that the sum over the stencil
    equals 1, ensuring strict conservation.

Usage
    \verbatim
    sourceSpreading
    {
        type    gaussian;
        sigma   0.01;   // Gaussian spreading width [m]
        nLayers 3;      // Topological stencil extent
    }
    \endverbatim

Note
    On processor boundaries only processor-local faces are included in
    the stencil, which truncates the Gaussian kernel.

SourceFiles
    gaussianSpreading.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_gaussianSpreading_H
#define regionFaModels_gaussianSpreading_H

#include "sourceSpreadingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                     Class gaussianSpreading Declaration
\*---------------------------------------------------------------------------*/

class gaussianSpreading
:
    public sourceSpreadingModel
{
    // Private Data

        //- Gaussian spreading width [m]
        scalar sigma_;

        //- Number of topological layers for BFS stencil
        label nLayers_;


    // Private Member Functions

        //- Templated implementation of spread
        template<class Type>
        void spreadImpl
        (
            Field<Type>& field,
            const scalarField& diameter
        ) const;

        //- No copy construct
        gaussianSpreading(const gaussianSpreading&) = delete;

        //- No copy assignment
        void operator=(const gaussianSpreading&) = delete;


public:

    //- Runtime type information
    TypeName("gaussian");


    // Constructors

        //- Construct from mesh and dictionary
        gaussianSpreading(const faMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~gaussianSpreading() = default;


    // Member Functions

        //- Spread scalar source field using Gaussian kernel
        virtual void spread
        (
            scalarField& field,
            const scalarField& diameter
        ) const;

        //- Spread vector source field using Gaussian kernel
        virtual void spread
        (
            vectorField& field,
            const scalarField& diameter
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
