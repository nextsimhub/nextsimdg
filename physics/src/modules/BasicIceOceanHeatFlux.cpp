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
    bool doPrint = (i == ModelArray::indexFromLocation(ModelArray::Type::H, {79,67}));
//    if (cice[i] > 0.) {
        if (doPrint) std::cerr << "BIOHF:" << " tf=" << tf[i] << " sst=" << sst[i] << " cpml=" << mlBulkCp[i] << std::endl;
        qio[i] = doOne(tf[i], sst[i], mlBulkCp[i], tst.step.seconds());
//    }
}

} /* namespace Nextsim */
