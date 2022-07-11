/*!
 * @file DummyOceanState.cpp
 *
 * @date May 10, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/DummyOceanState.hpp"

namespace Nextsim {
void DummyOceanState::setData(const ModelState& ms)
{
    sst = ms.at("sst");
    sss = ms.at("sss");
    mld = 10;
}
} /* namespace Nextsim */
