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

Class
    Foam::regionModels::areaSurfaceFilmModels::filmThicknessDiameter

Description
    Physical droplet-diameter model with \f$d=h\f$.  It is the baseline for
    direct droplets and drips and retains the legacy direct-transfer scaling
    when no explicit morphology model is selected.

Usage
    \verbatim
    diameterModels
    {
        <morphology>
        {
            // Mandatory entries
            type  filmThickness;
        }
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Reqd | Deflt
      type         | Type name: filmThickness           | word | yes  | -
    \endtable

SourceFiles
    filmThicknessDiameter.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_filmThicknessDiameter_H
#define regionFaModels_filmThicknessDiameter_H

#include "diameterModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                    Class filmThicknessDiameter Declaration
\*---------------------------------------------------------------------------*/

class filmThicknessDiameter
:
    public diameterModel
{
public:

    //- Runtime type information
    TypeName("filmThickness");


    // Constructors

        //- Construct from the film and model dictionary
        filmThicknessDiameter
        (
            const liquidFilmBase& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~filmThicknessDiameter() = default;


    // Member Functions

        //- Return the clipped local film thickness as a physical diameter
        virtual tmp<scalarField> diameter() const;

        //- Preserve legacy direct-transfer scaling
        virtual bool scalesWithSeparationRatio() const noexcept { return true; }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
