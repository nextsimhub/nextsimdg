/*!
 * @file BasicIceOceanHeatFlux.cpp
 *
 * @date Oct 19, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/BasicIceOceanHeatFlux.hpp"

namespace Nextsim {

static inline double doOne(double tBot, double sst, double mlBulkCp, double timeT)
{
    // Transfer rate depends on the mixed layer depth and the relaxation time scale
    return (sst - tBot) * mlBulkCp / timeT;
}

void BasicIceOceanHeatFlux::update(const TimestepTime& tst)
{
    overElements(std::bind(&BasicIceOceanHeatFlux::updateElement, this, std::placeholders::_1,
                     std::placeholders::_2),
        tst);
}

void BasicIceOceanHeatFlux::updateElement(size_t i, const TimestepTime& tst)
{
    // Use the timestep length as the relaxation time scale
    if (cice[i] > 0.) {
        qio[i] = doOne(tf[i], sst[i], mlBulkCp[i], tst.step.seconds());
    }
}

} /* namespace Nextsim */
