/*!
 * @file HiblerConcentration.cpp
 *
 * @date Nov 11, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/HiblerConcentration.hpp"

#include "include/NextsimPhysics.hpp"
#include "include/PhysicsData.hpp"
#include "include/PrognosticData.hpp"

namespace Nextsim {

double HiblerConcentration::h0 = 0;
double HiblerConcentration::phiM = 0.;

double HiblerConcentration::freeze(
    const PrognosticData& prog, PhysicsData& phys, NextsimPhysics& nsphys) const
{
    // Set the value of the reciprocal on the first invocation of freeze()
    static const double ooh0 = 1. / h0;
    return nsphys.newIce() * ooh0;
}

double HiblerConcentration::melt(
    const PrognosticData& prog, PhysicsData& phys, NextsimPhysics& nsphys) const
{
    if (prog.iceConcentration() >= 1) return 0;
    double del_hi = phys.updatedIceTrueThickness() - prog.iceTrueThickness();
    return del_hi * prog.iceConcentration() * phiM / prog.iceTrueThickness();
}

} /* namespace Nextsim */
