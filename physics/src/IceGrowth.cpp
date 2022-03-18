/*!
 * @file IceGrowth.cpp
 *
 * @date Mar 15, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/IceGrowth.hpp"
#include "include/VerticalIceGrowthModule.hpp"

namespace Nextsim {

const std::map<int, std::string> keyMap = {
    { IceGrowth::VERTICAL_GROWTH_KEY, "VerticalIceModel" },
    { IceGrowth::LATERAL_GROWTH_KEY, "LateralIceModel" },
};

void IceGrowth::configure()
{
    // Configure the vertical and lateral growth modules
    iVertical = std::move(Module::getInstance<IVerticalIceGrowth>());
}

void IceGrowth::update(const TimePoint& tsInitialTime)
{
    iVertical->update(tsInitialTime);
}
} /* namespace Nextsim */
