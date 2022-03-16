/*!
 * @file IceGrowth.cpp
 *
 * @date Mar 15, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/IceGrowth.hpp"

namespace Nextsim {

const std::map<int, std::string> keyMap = {
        {IceGrowth::VERTICAL_GROWTH_KEY, "VerticalIceModel"},
        {IceGrowth::LATERAL_GROWTH_KEY, "LateralIceModel"},
};

void IceGrowth::configure()
{

}
} /* namespace Nextsim */
