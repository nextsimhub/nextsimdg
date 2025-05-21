/*!
 * @file ERA5Atmosphere.cpp
 *
 * @date 21 May 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ERA5Atmosphere.hpp"

#include "include/Finalizer.hpp"
#include "include/NextsimModule.hpp"
#include "include/ParaGridIO.hpp"

#include <include/PhysicalBounds.hpp>

namespace Nextsim {

std::string ERA5Atmosphere::filePath;

bool ERA5Atmosphere::checkFieldsDefault = false;

std::string ERA5Atmosphere::pfx = "ERA5Atmosphere";
std::string ERA5Atmosphere::fileKey = pfx + ".file";
std::string ERA5Atmosphere::checkFieldsKey = pfx + ".check_fields";
std::string ERA5Atmosphere::fieldNamesKey = pfx + ".fields_names";
std::string ERA5Atmosphere::fieldNamesDefault
    = "tair,dew2m,pair,sw_in,lw_in,wind_speed,wind_u,wind_v";

static const std::map<int, std::string> keyMap = {
    { ERA5Atmosphere::FILEPATH_KEY, ERA5Atmosphere::fileKey },
    { ERA5Atmosphere::FIELDNAMES_KEY, ERA5Atmosphere::fieldNamesKey },
    { ERA5Atmosphere::CHECKFIELDS_KEY, ERA5Atmosphere::checkFieldsKey },
};

ERA5Atmosphere::ERA5Atmosphere()
    : fluxImpl(nullptr)
{
    getStore().registerArray(Protected::T_AIR, &tair, RO);
    getStore().registerArray(Protected::DEW_2M, &tdew, RO);
    getStore().registerArray(Protected::P_AIR, &pair, RO);
    getStore().registerArray(Protected::SW_IN, &sw_in, RO);
    getStore().registerArray(Protected::LW_IN, &lw_in, RO);
    getStore().registerArray(Protected::WIND_SPEED, &wind, RO);
}

ConfigurationHelp::HelpMap& ERA5Atmosphere::getHelpRecursive(HelpMap& map, bool getAll)
{
    auto options
        = getCheckingHelpList(fieldNamesKey, fieldNamesDefault, checkFieldsKey, checkFieldsDefault);
    options.push_back({ fileKey, ConfigType::STRING, {}, "", "",
        "Path to the processed NetCDF file providing the ERA5 forcings." });
    map[pfx] = options;
    Module::getHelpRecursive<IFluxCalculation>(map, getAll);

    return map;
}

void ERA5Atmosphere::configure()
{
    Finalizer::registerUnique(Module::finalize<IFluxCalculation>);

    filePath = Configured::getConfiguration(keyMap.at(FILEPATH_KEY), std::string());

    fluxImpl = &Module::getImplementation<IFluxCalculation>();
    tryConfigure(fluxImpl);

    if (Configured::getConfiguration(keyMap.at(CHECKFIELDS_KEY), checkFieldsDefault)) {
        setFieldsToCheck(fieldNamesDefault, pfx);
    }
}

ConfigMap ERA5Atmosphere::getConfiguration() const
{
    return {
        { keyMap.at(FILEPATH_KEY), filePath },
    };
}

void ERA5Atmosphere::update(const TimestepTime& tst)
{
    // TODO: Get more authoritative names for the forcings
    std::set<std::string> forcings
        = { "tair", "dew2m", "pair", "sw_in", "lw_in", "wind_speed", "u", "v" };

    ModelState state = ParaGridIO::readForcingTimeStatic(forcings, tst.start, filePath);
    tair = state.data.at("tair");
    tdew = state.data.at("dew2m");
    pair = state.data.at("pair");
    sw_in = state.data.at("sw_in");
    lw_in = state.data.at("lw_in");
    wind = state.data.at("wind_speed");
    uwind = state.data.at("u");
    vwind = state.data.at("v");
    snow = 0; // FIXME get snow data
    emp = 0; // FIXME get E - P data

    checkFields(tst);

    fluxImpl->update(tst);
}

void ERA5Atmosphere::setFilePath(const std::string& filePathIn) { filePath = filePathIn; }

void ERA5Atmosphere::setData(const ModelState::DataMap& ms)
{
    IAtmosphereBoundary::setData(ms);
    fluxImpl->setData(ms);
}

} /* namespace Nextsim */
