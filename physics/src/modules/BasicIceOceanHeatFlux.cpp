/*!
 * @file BasicIceOceanHeatFlux.cpp
 *
 * @date Oct 19, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/BasicIceOceanHeatFlux.hpp"

namespace Nextsim {

static inline double doOne(double tBot, double sst, double mlBulkCp, double ts)
{
    // Transfer rate depends on the mixed layer depth an d a timescale. Here, it is the timestep
    return (sst - tBot) * mlBulkCp / ts;
}

void BasicIceOceanHeatFlux::update(const TimestepTime& tst)
{
    overElements(std::bind(&BasicIceOceanHeatFlux::updateElement, this, std::placeholders::_1,
                     std::placeholders::_2),
        tst);
}

void BasicIceOceanHeatFlux::updateElement(size_t i, const TimestepTime& tst)
{
    qio[i] = doOne(tf[i], sst[i], mlBulkCp[i], tst.step.seconds());
}

} /* namespace Nextsim */
