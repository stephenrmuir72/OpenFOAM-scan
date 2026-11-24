/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2009 Frank Bos, TU Delft
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

#include "RBFFunction.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(RBFFunction);
    defineRunTimeSelectionTable(RBFFunction, dictionary);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::autoPtr<Foam::RBFFunction> Foam::RBFFunction::New
(
    const word& RBFname,
    const dictionary& dict
)
{
    Info<< "Radial Basis Function interpolation: " << RBFname << endl;

    auto* ctorPtr = dictionaryConstructorTable(RBFname);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "Unknown RBF function type",
            RBFname,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    const auto& coeffs =
        dict.optionalSubDict(RBFname + "Coeffs", keyType::LITERAL);

    return autoPtr<RBFFunction>(ctorPtr(coeffs));
}


// ************************************************************************* //
