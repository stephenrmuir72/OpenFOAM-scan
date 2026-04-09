DebugVar("LLLL");

CorrectPhi
(
    U,
    phi,
    p,
    dimensionedScalar("rAUf", dimTime, 1),
    geometricZeroField(),
    pimple
);

DebugVar("MMMM");

#include "continuityErrs.H"

DebugVar("NNNN");
