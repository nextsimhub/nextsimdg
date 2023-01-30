/*!
 * @file BasicIceOceanHeatFlux.cpp
 *
 * @date Oct 19, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/BasicIceOceanHeatFlux.hpp"

namespace Nextsim {

double BasicIceOceanHeatFlux::timeT = 2592000; // Relaxation time, s. Defaults to 30 days.

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
    qio[i] = doOne(tf[i], sst[i], mlBulkCp[i], timeT);
}

} /* namespace Nextsim */
