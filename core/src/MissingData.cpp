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
static const double missingDefault = -0x1p300;
double MissingData::m_value = missingDefault;
void MissingData::configure()
{
    m_value = Configured::getConfiguration(keyMap.at(MissingData::MISSINGVALUE_KEY), missingDefault);
}

MissingData::HelpMap& MissingData::getHelpText(HelpMap& map, bool getAll)
{
    map["MissingData"] = {
            { keyMap.at(MISSINGVALUE_KEY), ConfigType::NUMERIC, {"-∞", "∞"}, "-2³⁰⁰", "",
                    "Missing data indicator used for input and output."
            },
    };
    return map;
}
} /* namespace Nextsim */
