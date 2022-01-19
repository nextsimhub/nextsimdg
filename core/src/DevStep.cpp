/*!
 * @file DevStep.cpp
 *
 * @date Jan 12, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/DevStep.hpp"

namespace Nextsim {

void DevStep::iterate(const Iterator::Duration& dt)
{
    PrognosticData::setTimestep(dt);
    for (pStructure->cursor = 0; pStructure->cursor; ++pStructure->cursor) {
        ElementData data = *pStructure->cursor;
        data.updateDerivedData(data, data, data);
        data.calculate(data, data, data);
    }
}

} /* namespace Nextsim */
