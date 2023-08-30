/*!
 * @file ERA5Atmosphere.cpp
 *
 * @date Nov 25, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ERA5Atmosphere.hpp"

#include "include/Module.hpp"
#include "include/ParaGridIO.hpp"

namespace Nextsim {

std::string ERA5Atmosphere::filePath;

static const std::string pfx = "ERA5Atmosphere";
static const std::string fileKey = pfx + ".file";

template <>
const std::map<int, std::string> Configured<ERA5Atmosphere>::keyMap = {
    { ERA5Atmosphere::FILEPATH_KEY, fileKey },
};

ERA5Atmosphere::ERA5Atmosphere()
    : fluxImpl(0)
{
    getStore().registerArray(Protected::T_AIR, &tair);
    getStore().registerArray(Protected::DEW_2M, &tdew);
    getStore().registerArray(Protected::P_AIR, &pair);
    getStore().registerArray(Protected::SW_IN, &sw_in);
    getStore().registerArray(Protected::LW_IN, &lw_in);
    getStore().registerArray(Protected::WIND_SPEED, &wind);
}

ConfigurationHelp::HelpMap& ERA5Atmosphere::getHelpRecursive(HelpMap& map, bool getAll)
{
    map[pfx] = {
        { fileKey, ConfigType::STRING, {}, "", "",
            "Path to the processed NetCDF file providing the ERA5 forcings." },
    };
    Module::getHelpRecursive<IFluxCalculation>(map, getAll);

    return map;
}

void ERA5Atmosphere::configure()
{
    filePath = Configured::getConfiguration(keyMap.at(FILEPATH_KEY), std::string());

    fluxImpl = &Module::getImplementation<IFluxCalculation>();
    tryConfigure(fluxImpl);
}

void ERA5Atmosphere::update(const TimestepTime& tst)
{
    // TODO: Get more authoritative names for the forcings
    std::set<std::string> forcings = { "tair", "dew2m", "pair", "sw_in", "lw_in", "wind_speed", "u", "v" };

    ModelState state = ParaGridIO::readForcingTimeStatic(forcings, tst.start, filePath);
    tair = state.data.at("tair");
    tdew = state.data.at("dew2m");
    pair = state.data.at("pair");
    sw_in = state.data.at("sw_in");
    lw_in = state.data.at("lw_in");
    wind = state.data.at("wind_speed");
    uwind = state.data.at("u");
    vwind = state.data.at("v");

    fluxImpl->update(tst);
}

void ERA5Atmosphere::setFilePath(const std::string& filePathIn) { filePath = filePathIn; }

void ERA5Atmosphere::setData(const ModelState::DataMap& ms)
{
    IAtmosphereBoundary::setData(ms);
    fluxImpl->setData(ms);
}

} /* namespace Nextsim */
