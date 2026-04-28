/*---------------------------------------------------------------------------*\

License
    This file is part of BRNSFoam, derived from GeoChemFoam.

    BRNSFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version. See <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "makeTwoPhaseMultiComponentTransportMixture.H"

#include "basicTwoPhaseMultiComponentTransportMixture.H"
#include "twoPhaseMultiComponentTransportMixture.H"

#include "inertMultiComponentMixture.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeTwoPhaseMultiComponentTransportMixture
(
    twoPhaseMultiComponentTransportMixture,
    inertMultiComponentMixture,
    inertMultiComponentMixture
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
