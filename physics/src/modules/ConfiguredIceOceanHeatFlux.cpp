/*!
 * @file ConfiguredIceOceanHeatFlux.cpp
 *
 * @date 5 Jul 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ConfiguredIceOceanHeatFlux.hpp"

#include <string>

#include <iostream>

namespace Nextsim {

double ConfiguredIceOceanHeatFlux::qio0 = 0;

static const std::string pfx = "ConfiguredIceOceanHeatFlux";
static const std::string qioKey = pfx + ".qio";

template <>
const std::map<int, std::string> Configured<ConfiguredIceOceanHeatFlux>::keyMap = {
    { ConfiguredIceOceanHeatFlux::QIO_KEY, qioKey },
};

ConfiguredIceOceanHeatFlux::ConfiguredIceOceanHeatFlux()
    : IIceOceanHeatFlux()
    , qioConf(qio0)
{
}

ConfigurationHelp::HelpMap& ConfiguredIceOceanHeatFlux::getHelpRecursive(HelpMap& map, bool getAll)
{
    map[pfx] = {
        { qioKey, ConfigType::NUMERIC, { "-∞", "∞" }, std::to_string(qio0), "W m⁻²",
            "Ice-ocean heat flux (+ve upward)." },
    };
    return map;
}

void ConfiguredIceOceanHeatFlux::configure()
{
    qioConf = Configured<ConfiguredIceOceanHeatFlux>::getConfiguration(
        Configured<ConfiguredIceOceanHeatFlux>::keyMap.at(QIO_KEY), qio0);

    std::cerr << "Configuring ice ocean heat flux to " << qioConf << std::endl;
}

void ConfiguredIceOceanHeatFlux::setData(const ModelState::DataMap&)
{
    qio.data() = qioConf;
    std::cerr << "Setting ice ocean heat flux to " << qioConf << std::endl;
}
} /* namespace Nextsim */
