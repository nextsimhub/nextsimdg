/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ERA5Atmosphere.hpp"

#include "include/Finalizer.hpp"
#include "include/NextsimModule.hpp"
#include "include/ParaGridIO.hpp"

namespace Nextsim {

std::string ERA5Atmosphere::filePath;

static const std::string pfx = "ERA5Atmosphere";
static const std::string fileKey = pfx + ".file";

static const std::map<int, std::string> keyMap = {
    { ERA5Atmosphere::FILEPATH_KEY, fileKey },
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
    map[pfx] = {
        { fileKey, ConfigType::STRING, {}, "", "",
            "Path to the processed NetCDF file providing the ERA5 forcings." },
    };
    Module::getHelpRecursive<IFluxCalculation>(map, getAll);

    return map;
}

void ERA5Atmosphere::configure()
{
    Finalizer::registerUnique(Module::finalize<IFluxCalculation>);

    filePath = Configured::getConfiguration(keyMap.at(FILEPATH_KEY), std::string());

    fluxImpl = &Module::getImplementation<IFluxCalculation>();
    tryConfigure(fluxImpl);

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

    // Read ERA5 forcings at the top of the hour
    if (tst.start.format(TimePoint::msFormat) == "T00:00Z") {
        forcingState = ParaGridIO::readForcingTimeStatic(forcings, tst.start, filePath);
    }
    tair = forcingState.data.at("tair");
    tdew = forcingState.data.at("dew2m");
    pair = forcingState.data.at("pair");
    sw_in = forcingState.data.at("sw_in");
    lw_in = forcingState.data.at("lw_in");
    wind = forcingState.data.at("wind_speed");
    uwind = forcingState.data.at("u");
    vwind = forcingState.data.at("v");
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
