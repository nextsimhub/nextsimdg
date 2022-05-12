/*!
 * @file DevStep.cpp
 *
 * @date Jan 12, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/DevStep.hpp"

namespace Nextsim {

void DevStep::iterate(const TimestepTime& tst)
{
    pData->update(tst);
}

} /* namespace Nextsim */
