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

\*---------------------------------------------------------------------------*/

#include "morphologyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::regionModels::areaSurfaceFilmModels::morphologyModel>
Foam::regionModels::areaSurfaceFilmModels::morphologyModel::New
(
    const liquidFilmBase& film,
    const dictionary& dict
)
{
    const word modelType(dict.get<word>("type"));

    Info<< "    Selecting film-separation morphology model "
        << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "morphologyModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<morphologyModel>(ctorPtr(film, dict));
}


Foam::autoPtr<Foam::regionModels::areaSurfaceFilmModels::morphologyModel>
Foam::regionModels::areaSurfaceFilmModels::morphologyModel::NewIfPresent
(
    const liquidFilmBase& film,
    const dictionary& filmDict
)
{
    if (!filmDict.found("morphology"))
    {
        return nullptr;
    }

    return New(film, filmDict.subDict("morphology"));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
