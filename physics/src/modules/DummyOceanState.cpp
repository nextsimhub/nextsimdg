/*!
 * @file DummyOceanState.cpp
 *
 * @date May 10, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/DummyOceanState.hpp"

namespace Nextsim {
void DummyOceanState::setData(const ModelState::DataMap& ms)
{
    sst = -1.5; // Liquid ocean, solid ice
    sss = 32; // Very rough average for the Arctic Ocean
    mld = 10;
}
} /* namespace Nextsim */
