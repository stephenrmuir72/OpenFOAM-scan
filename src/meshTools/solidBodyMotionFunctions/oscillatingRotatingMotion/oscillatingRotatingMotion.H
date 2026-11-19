/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::solidBodyMotionFunctions::oscillatingRotatingMotion

Description
    SolidBodyMotionFvMesh 6DoF motion function. Oscillating rotation.

    The governing equation for the rotation is given by:
    \f[
        r = A*\sin(B*(t + C)) + D
    \f]

    where:
    \vartable
        r    | Rotation angles [deg] (vector)
        A    | Rotation angles [deg] (vector)
        B    | Radial velocity [rad/s] (scalar)
        C    | Phase-shift to left [s] (scalar)
        D    | Angle shift [deg] (vector)
    \endvartable

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    oscillation
    {
        // Mandatory entries
        solidBodyMotionFunction oscillatingRotatingMotion;

        oscillatingRotatingMotionCoeffs
        {
            // Mandatory entries
            origin          <vector>;
            amplitude       <Function1<vector>>;
            omega           <Function1<scalar>>;

            // Optional entries
            phaseShift      <Function1<scalar>>;
            angleShift      <Function1<vector>>;
        }
    }

    where the entries mean:
    \table
      Property  | Description          | Type                | Reqd | Deflt
      solidBodyMotionFunction          | Type name: oscillatingRotatingMotion <!--
                -->                    | word                | yes  | -
      origin    | Centre of rotation   | <vector>            | yes  | -
      amplitude | Rotation angles      | Function1\<vector\> | yes  | -
      omega     | Radial velocity      | Function1\<scalar\> | yes  | -
      phaseShift | Phase-shift to left | Function1\<scalar\> | no   | -
      angleShift | Angle offset        | Function1\<vector\> | no   | -
    \endtable

SourceFiles
    oscillatingRotatingMotion.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_solidBodyMotionFunctions_oscillatingRotatingMotion_H
#define Foam_solidBodyMotionFunctions_oscillatingRotatingMotion_H

#include "solidBodyMotionFunction.H"
#include "primitiveFields.H"
#include "point.H"
#include "Function1.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                  Class oscillatingRotatingMotion Declaration
\*---------------------------------------------------------------------------*/

class oscillatingRotatingMotion
:
    public solidBodyMotionFunction
{
    // Private Data

        //- Centre of rotation
        point origin_;

        //- Radial velocity [rad/s]
        autoPtr<Function1<scalar>> omegaPtr_;

        //- Amplitude == rotation angles [deg]
        autoPtr<Function1<vector>> amplitudePtr_;

        //- Phase shift [s] (optional)
        autoPtr<Function1<scalar>> phaseShiftPtr_;

        //- Angle shift [deg] (optional)
        autoPtr<Function1<vector>> angleShiftPtr_;


public:

    //- Runtime type information
    TypeName("oscillatingRotatingMotion");


    // Generated Methods

        //- No copy construct
        oscillatingRotatingMotion(const oscillatingRotatingMotion&) = delete;

        //- No copy assignment
        void operator=(const oscillatingRotatingMotion&) = delete;


    // Constructors

        //- Construct from components
        oscillatingRotatingMotion
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );

        //- Construct and return a clone
        virtual autoPtr<solidBodyMotionFunction> clone() const
        {
            return autoPtr<solidBodyMotionFunction>
            (
                new oscillatingRotatingMotion
                (
                    SBMFCoeffs_,
                    time_
                )
            );
        }


    //- Destructor
    virtual ~oscillatingRotatingMotion() = default;


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
