/*!
 * @file DevStep.cpp
 *
 * @date Jan 12, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/DevStep.hpp"
#include "include/IPrognosticUpdater.hpp"
#include "include/PrognosticElementData.hpp"

namespace Nextsim {

void DevStep::iterate(const TimestepTime& tst)
{
//    PrognosticElementData::setTimestep(tst.step);
//    for (pStructure->cursor = 0; pStructure->cursor; ++pStructure->cursor) {
//        auto& data = *pStructure->cursor;
//        data.updateDerivedData(data, data, data);
//        data.calculate(data, data, data);
//        data.updateAndIntegrate(data);
//    }
    pData->update(tst);
}

} /* namespace Nextsim */
