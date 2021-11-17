/*!
 * @file HiblerConcentration.cpp
 *
 * @date Nov 11, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/HiblerConcentration.hpp"

#include "include/ModuleLoader.hpp"
#include "include/NextsimPhysics.hpp"
#include "include/PhysicsData.hpp"
#include "include/PrognosticData.hpp"

namespace Nextsim {

double HiblerConcentration::h0 = 0;
double HiblerConcentration::phiM = 0.;

static std::map<int, std::string> Configured<HiblerConcentration>::keyMap = {
    { HiblerConcentration::H0_KEY, "Hibler.h0" },
    { HiblerConcentration::PHIM_KEY, "Hibler.phiM" },
};

void HiblerConcentration::configure()
{
    ModuleLoader& loader = ModuleLoader::getLoader();

    h0 = Configured::getConfiguration(keyMap.at(H0_KEY), 0.25);
    phiM = Configured::getConfiguration(keyMap.at(PHIM_KEY), 0.5);
}

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
    if (prog.iceConcentration() >= 1)
        return 0;
    double del_hi = phys.updatedIceTrueThickness() - prog.iceTrueThickness();
    return del_hi * prog.iceConcentration() * phiM / prog.iceTrueThickness();
}

} /* namespace Nextsim */
