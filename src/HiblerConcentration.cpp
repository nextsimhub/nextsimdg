/*!
 * @file HiblerConcentration.cpp
 *
 * @date Nov 11, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/HiblerConcentration.hpp"

#include "include/PrognosticData.hpp"
#include "include/NextsimPhysics.hpp"

namespace Nextsim {

double HiblerConcentration::h0 = 0;

double HiblerConcentration::freeze(const PrognosticData&, NextsimPhysics&) const
{

    return 0;
}

double HiblerConcentration::melt(const PrognosticData&, NextsimPhysics&) const
{
    return 0;
}

} /* namespace Nextsim */
