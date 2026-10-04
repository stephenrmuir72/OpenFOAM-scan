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

#include "temporalSpreading.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::regionModels::areaSurfaceFilmModels::temporalSpreading>
Foam::regionModels::areaSurfaceFilmModels::temporalSpreading::New
(
    const faMesh& mesh,
    const dictionary& dict
)
{
    const word modelType(dict.get<word>("type"));

    Info<< "    Selecting temporalSpreading " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "temporalSpreading",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<temporalSpreading>(ctorPtr(mesh, dict));
}


Foam::autoPtr<Foam::regionModels::areaSurfaceFilmModels::temporalSpreading>
Foam::regionModels::areaSurfaceFilmModels::temporalSpreading::NewIfPresent
(
    const faMesh& mesh,
    const dictionary& filmDict
)
{
    if (!filmDict.found("temporalSpreading"))
    {
        localIOdictionary state(stateIO(mesh), word::null);

        if
        (
            state.found("stateVersion")
         && !state.found("sparseFaces")
        )
        {
            FatalErrorInFunction
                << "Incomplete temporalSpreadingState: missing sparseFaces"
                << abort(FatalError);
        }

        if
        (
            state.found("stateVersion")
         &&
            (
                state.get<label>("stateVersion") == 3
             || state.get<label>("stateVersion") == 4
            )
         && !state.found("sparseFaceLabels")
        )
        {
            FatalErrorInFunction
                << "Incomplete temporalSpreadingState: missing "
                << "sparseFaceLabels"
                << abort(FatalError);
        }

        if (state.found("eventIds") || state.found("sparseFaces"))
        {
            const labelList eventIds
            (
                state.found("eventIds")
              ? state.get<labelList>("eventIds") : labelList()
            );
            const labelList sparseFaces
            (
                state.found("sparseFaces")
              ? state.get<labelList>("sparseFaces") : labelList()
            );

            if (!eventIds.empty() || !sparseFaces.empty())
            {
                FatalErrorInFunction
                    << "Cannot remove temporalSpreading while "
                    << eventIds.size()
                    << " dense event(s) and " << sparseFaces.size()
                    << " sparse contribution(s) remain in "
                    << "temporalSpreadingState. Configure the saved model "
                    << "until the events have released, or remove the "
                    << "pending state deliberately."
                    << abort(FatalError);
            }
        }

        return nullptr;
    }

    return New(mesh, filmDict.subDict("temporalSpreading"));
}


// ************************************************************************* //
