/*!
 * @file ERA5Atmosphere.cpp
 *
 * @date 30 May 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ERA5Atmosphere.hpp"

#include "include/Finalizer.hpp"
#include "include/NextsimModule.hpp"
#include "include/ParaGridIO.hpp"

namespace Nextsim {

std::string ERA5Atmosphere::filePath;

bool ERA5Atmosphere::checkFieldsDefault = false;

std::string ERA5Atmosphere::pfx = "ERA5Atmosphere";
std::string ERA5Atmosphere::fileKey = pfx + ".file";
std::string ERA5Atmosphere::checkFieldsKey = pfx + ".check_fields";
std::string ERA5Atmosphere::fieldNamesKey = pfx + ".fields_names";

static const std::map<int, std::string> keyMap = {
    { ERA5Atmosphere::FILEPATH_KEY, ERA5Atmosphere::fileKey },
    { ERA5Atmosphere::CHECKFIELDS_KEY, ERA5Atmosphere::checkFieldsKey },
};

ERA5Atmosphere::ERA5Atmosphere()
    : fluxImpl(nullptr)
    , tair(ModelArray::Type::H, { -100, 100 })
    , tdew(ModelArray::Type::H, { -100, 100 })
    , pair(ModelArray::Type::H, { 500e2, 2000e2 })
    , sw_in(ModelArray::Type::H, { -1e-6, 1e4 })
    , lw_in(ModelArray::Type::H, { -1e-6, 1e4 })
    , wind(ModelArray::Type::H, { 0, 100 })
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
    map[pfx] = { { fileKey, ConfigType::STRING, {}, "", "",
                     "Path to the processed NetCDF file providing the ERA5 forcings." },
        { checkFieldsKey, ConfigType::BOOLEAN, { "true", "false" },
            ConfigurationHelp::toString(checkFieldsDefault), "",
            "Set to true to check if the main module variables fall within a reasonable physical "
            "range." } };

    Module::getHelpRecursive<IFluxCalculation>(map, getAll);

    return map;
}

void ERA5Atmosphere::configure()
{
    Finalizer::registerUnique(Module::finalize<IFluxCalculation>);

    filePath = Configured::getConfiguration(keyMap.at(FILEPATH_KEY), std::string());

    fluxImpl = &Module::getImplementation<IFluxCalculation>();
    tryConfigure(fluxImpl);

    boolCheckFields = Configured::getConfiguration(keyMap.at(CHECKFIELDS_KEY), checkFieldsDefault);

    addChecks({
        { "tair", &tair },
        { "tdew", &tdew },
        { "pair", &pair },
        { "sw_in", &sw_in },
        { "lw_in", &lw_in },
        { "wind", &wind },
    });
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
    rain = 0; // FIXME get rain data

    fluxImpl->update(tst);

    try {
        checkFields();
    } catch (const std::exception& e) {
        throw std::runtime_error("ERA5Atmosphere:update: " + std::string(e.what()));
    }
}

void ERA5Atmosphere::setFilePath(const std::string& filePathIn) { filePath = filePathIn; }

void ERA5Atmosphere::setData(const ModelState::DataMap& ms)
{
    IAtmosphereBoundary::setData(ms);
    fluxImpl->setData(ms);
}

} /* namespace Nextsim */
