/*!
 * @file ModelData.cpp
 *
 * @date Jul 20, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ModelData.hpp"

namespace Nextsim {

void ModelData::updateAtmosOceanState(const TimestepTime& tst)
{
    aoState.update(tst);
}

void ModelData::stepPrognosticData(const TimestepTime& tst)
{
    pData.update(tst);

    metadata.incrementTime(tst.step);
}

void ModelData::exportData(const TimestepTime& tst) const
{
    // Nothing yet
}

ModelState ModelData::getModelState() const
{
    return pData.getStateRecursive(true);
}

} /* namespace Nextsim */
