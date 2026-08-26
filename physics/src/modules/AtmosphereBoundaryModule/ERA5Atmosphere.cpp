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
    , tairAccessor(getStore(), RO, ModelArray::Type::H, std::pair(-100.0, 100.0))
    , tdewAccessor(getStore(), RO, ModelArray::Type::H, std::pair(-100.0, 100.0))
    , pairAccessor(getStore(), RO, ModelArray::Type::H, std::pair(500e2, 2000e2))
    , sw_inAccessor(getStore(), RO, ModelArray::Type::H, std::pair(-1e-3, 1e4))
    , lw_inAccessor(getStore(), RO, ModelArray::Type::H, std::pair(-1e-3, 1e4))
    , windAccessor(getStore(), RO, ModelArray::Type::H, std::pair(0.0, 100.0))
{
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
        { "tair", tairAccessor },
        { "tdew", tdewAccessor },
        { "pair", pairAccessor },
        { "sw_in", sw_inAccessor },
        { "lw_in", lw_inAccessor },
        { "wind", windAccessor },
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
    if (std::fmod((tst.start - TimePoint()).seconds(), 3600.) == 0.0_ft) {
        forcingState = ParaGridIO::readForcingTimeStatic(forcings, tst.start, filePath);
    }
    tairAccessor.getHostRW() = forcingState.data.at("tair");
    tdewAccessor.getHostRW() = forcingState.data.at("dew2m");
    pairAccessor.getHostRW() = forcingState.data.at("pair");
    sw_inAccessor.getHostRW() = forcingState.data.at("sw_in");
    lw_inAccessor.getHostRW() = forcingState.data.at("lw_in");
    windAccessor.getHostRW() = forcingState.data.at("wind_speed");
    uwindAccessor.getHostRW() = forcingState.data.at("u");
    vwindAccessor.getHostRW() = forcingState.data.at("v");
    snowAccessor.getHostRW() = 0; // FIXME get snow data
    rainAccessor.getHostRW() = 0; // FIXME get rain data

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
