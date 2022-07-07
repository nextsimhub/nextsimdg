/*!
 * @file MissingData.cpp
 *
 * @date Jun 14, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/MissingData.hpp"

namespace Nextsim {

template <>
const std::map<int, std::string> Configured<MissingData>::keyMap = {
    { MissingData::MISSINGVALUE_KEY, "model.missing_value" },
};
double MissingData::m_value = -0x1p300;
void MissingData::configure()
{
    m_value = Configured::getConfiguration(keyMap.at(MissingData::MISSINGVALUE_KEY), m_value);
}
} /* namespace Nextsim */
