/*!
 * @file HiblerSpread.cpp
 *
 * @date Apr 5, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/HiblerSpread.hpp"

namespace Nextsim {

double HiblerSpread::h0 = 0;
double HiblerSpread::phiM = 0;

template <>
const std::map<int, std::string> Configured<HiblerSpread>::keyMap = {
    { HiblerSpread::H0_KEY, "Hibler.h0" },
    { HiblerSpread::PHIM_KEY, "Hibler.phiM" },
};

void HiblerSpread::configure()
{
    h0 = Configured::getConfiguration(keyMap.at(H0_KEY), 0.25);
    phiM = Configured::getConfiguration(keyMap.at(PHIM_KEY), 0.5);
}

void HiblerSpread::freeze(const TimestepTime& tstep, double hice, double hsnow, double deltaHi,
    double newIce, double& cice, double& qow, double& deltaCfreeze)
{
    static const double ooh0 = 1. / h0;
    deltaCfreeze = newIce * ooh0;
}

void HiblerSpread::melt(const TimestepTime& tstep, double hice, double hsnow, double deltaHi,
    double& cice, double& qow, double& deltaCmelt)
{
    if (cice < 1) {
        deltaCmelt = deltaHi * cice * phiM / hice;
    }
}

}
