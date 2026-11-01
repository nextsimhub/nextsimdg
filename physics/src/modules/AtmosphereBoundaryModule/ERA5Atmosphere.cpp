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

    fluxImpl = Module::getInstance<IFluxCalculation>();
    tryConfigure(*fluxImpl);

    forcingState = Module::getInstance<IInput>();

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
    forcingState->update(tst.start);

    // Fetch the raw data fields
    tairAccessor.getHostRW() = forcingState->getField(tAirName);
    tdewAccessor.getHostRW() = forcingState->getField(dew2mName);
    pairAccessor.getHostRW() = forcingState->getField(pAirName);
    sw_inAccessor.getHostRW() = forcingState->getField(swInName);
    lw_inAccessor.getHostRW() = forcingState->getField(lwInName);
    uwindAccessor.getHostRW() = forcingState->getField(uName);
    vwindAccessor.getHostRW() = forcingState->getField(vName);

    // These ifs and warnings are for backwards compatability
    if (forcingState->isValid(snowName)) {
        snowAccessor.getHostRW() = forcingState->getField(snowName);
    } else {
        snowAccessor.getAutoRW() = 0._ft;
        // Only warn once
        if (warn.count(snowName) == 0) {
            Logged::warning("ERA5Atmosphere::update: No snow fall in forcing file, but I continue "
                            "anyway with zero snow fall.\n");
            warn[snowName] = true;
        }
    }

    if (forcingState->isValid(rainName)) {
        rainAccessor.getHostRW() = forcingState->getField(rainName);
    } else {
        rainAccessor.getAutoRW() = 0._ft;
        // Only warn once
        if (warn.count(rainName) == 0) {
            Logged::warning("ERA5Atmosphere::update: No rain fall in forcing file, but I continue "
                            "anyway with zero rain fall.\n");
            warn[rainName] = true;
        }
    }

    // Postprocess the data fields
    auto& tair = tairAccessor.getAutoRW();
    auto& tdew = tdewAccessor.getAutoRW();
    auto& rain = rainAccessor.getAutoRW();
    const auto& snow = snowAccessor.getAutoRO();
    auto& wind = windAccessor.getAutoRW();
    const auto& uWind = uwindAccessor.getAutoRO();
    const auto& vWind = vwindAccessor.getAutoRO();
    overElementsAuto(OVER_ELEMENTS_LAMBDA(const ElementIndex i) {
        // Temperatures in Celsius
        tair[i] -= Water::Tf;
        tdew[i] -= Water::Tf;

        // Rain from ERA5 is actually total precipitation - so subtract the snowfall
        rain[i] -= snow[i];

        // Calculate wind speed from components
        wind[i] = std::sqrt(uWind[i] * uWind[i] + vWind[i] * vWind[i]);
    });

    fluxImpl->update(tst);

    try {
        checkFields();
    } catch (const std::exception& e) {
        throw std::runtime_error("ERA5Atmosphere::update:: " + std::string(e.what()));
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
    forcingState->initialise(metadata.startTime(), filePath, ncLonName, ncLatName, ncTimeName,
        forcings, vectors, state.data.at(longitudeName), state.data.at(latitudeName));
}

} /* namespace Nextsim */
