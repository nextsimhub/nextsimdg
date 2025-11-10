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
    , sw_inAccessor(getStore(), RO, ModelArray::Type::H, std::pair(-1e-6, 1e4))
    , lw_inAccessor(getStore(), RO, ModelArray::Type::H, std::pair(-1e-6, 1e4))
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
        { "tair", &tairAccessor.getHostRO() },
        { "tdew", &tdewAccessor.getHostRO() },
        { "pair", &pairAccessor.getHostRO() },
        { "sw_in", &sw_inAccessor.getHostRO() },
        { "lw_in", &lw_inAccessor.getHostRO() },
        { "wind", &windAccessor.getHostRO() },
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
    tairAccessor.getHostRW() = state.data.at("tair");
    tdewAccessor.getHostRW() = state.data.at("dew2m");
    pairAccessor.getHostRW() = state.data.at("pair");
    sw_inAccessor.getHostRW() = state.data.at("sw_in");
    lw_inAccessor.getHostRW() = state.data.at("lw_in");
    windAccessor.getHostRW() = state.data.at("wind_speed");
    uwindAccessor.getHostRW() = state.data.at("u");
    vwindAccessor.getHostRW() = state.data.at("v");
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
