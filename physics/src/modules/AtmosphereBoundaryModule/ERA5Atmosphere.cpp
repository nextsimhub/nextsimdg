/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ERA5Atmosphere.hpp"

#include "include/Finalizer.hpp"
#include "include/NextsimModule.hpp"
#include "include/ParaGridIO.hpp"
#include "include/constants.hpp"

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

    fluxImpl = std::move(Module::getInstance<IFluxCalculation>());
    tryConfigure(*fluxImpl);

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
    forcingState.update(tst.start);

    tairAccessor.getHostRW() = forcingState.getField(tAirName) - Water::Tf;
    tdewAccessor.getHostRW() = forcingState.getField(dew2mName) - Water::Tf;
    pairAccessor.getHostRW() = forcingState.getField(pAirName);
    sw_inAccessor.getHostRW() = forcingState.getField(swInName);
    lw_inAccessor.getHostRW() = forcingState.getField(lwInName);
    uwindAccessor.getHostRW() = forcingState.getField(uName);
    vwindAccessor.getHostRW() = forcingState.getField(vName);
    snowAccessor.getHostRW() = forcingState.getField(snowName);
    rainAccessor.getHostRW() = forcingState.getField(rainName);

    windAccessor.getHostRW()
        = (uwindAccessor.getHostRW().data().pow(2) + vwindAccessor.getHostRW().data().pow(2))
              .sqrt();

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

    ModelState state;
    const ModelMetadata& metadata = ModelMetadata::getInstance();
    metadata.affixCoordinates(state);
    forcingState.setData(metadata.startTime(), filePath, ncLonName, ncLatName, ncTimeName, forcings,
        vectors, state.data[longitudeName], state.data[latitudeName]);
}

} /* namespace Nextsim */
