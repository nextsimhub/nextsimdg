/*!
 * @author  Einar Ólason <einar.olason@nersc.no>
 */

#include "include/CARRA2Atmosphere.hpp"

#include "FluxCalculationModule/include/FiniteElementSpecHum.hpp"
#include "include/Finalizer.hpp"
#include "include/NextsimModule.hpp"
#include "include/ParaGridIO.hpp"
#include "include/constants.hpp"

namespace Nextsim {

std::string CARRA2Atmosphere::filePath;

static const std::string pfx = "CARRA2Atmosphere";
static const std::string fileKey = pfx + ".file";

static const std::map<int, std::string> keyMap = {
    { CARRA2Atmosphere::FILEPATH_KEY, fileKey },
};

CARRA2Atmosphere::CARRA2Atmosphere()
    : tairAccessor(getStore(), RO, ModelArray::Type::H, std::pair(-100.0, 100.0))
    , spHumAccessor(getStore(), RO, ModelArray::Type::H, std::pair(0.0, 1.0))
    , pairAccessor(getStore(), RO, ModelArray::Type::H, std::pair(500e2, 2000e2))
    , sw_inAccessor(getStore(), RO, ModelArray::Type::H, std::pair(-1e-3, 1e4))
    , lw_inAccessor(getStore(), RO, ModelArray::Type::H, std::pair(-1e-3, 1e4))
    , windAccessor(getStore(), RO, ModelArray::Type::H, std::pair(0.0, 100.0))
    , fluxImpl(nullptr)
{
}

ConfigurationHelp::HelpMap& CARRA2Atmosphere::getHelpRecursive(HelpMap& map, bool getAll)
{
    map[pfx] = {
        { fileKey, ConfigType::STRING, {}, "", "",
            "Path to the processed NetCDF file providing the CARRA2 forcings." },
    };
    Module::getHelpRecursive<IFluxCalculation>(map, getAll);

    return map;
}

void CARRA2Atmosphere::configure()
{
    Finalizer::registerUnique(Module::finalize<IFluxCalculation>);

    filePath = Configured::getConfiguration(keyMap.at(FILEPATH_KEY), std::string());

    fluxImpl = Module::getInstance<IFluxCalculation>();
    tryConfigure(*fluxImpl);

    addChecks({
        { "tair", tairAccessor },
        { "pair", pairAccessor },
        { "sw_in", sw_inAccessor },
        { "lw_in", lw_inAccessor },
        { "wind", windAccessor },
    });
}

ConfigMap CARRA2Atmosphere::getConfiguration() const
{
    return {
        { keyMap.at(FILEPATH_KEY), filePath },
    };
}

void CARRA2Atmosphere::update(const TimestepTime& tst)
{
    forcingState->update(tst.start);

    // Fetch the raw data fields
    tairAccessor.getHostRW() = forcingState->getField(tAirName);
    spHumAccessor.getHostRW() = forcingState->getField(sphumName);
    pairAccessor.getHostRW() = forcingState->getField(pAirName);
    sw_inAccessor.getHostRW() = forcingState->getField(swInName);
    lw_inAccessor.getHostRW() = forcingState->getField(lwInName);
    uwindAccessor.getHostRW() = forcingState->getField(uName);
    vwindAccessor.getHostRW() = forcingState->getField(vName);
    snowAccessor.getHostRW() = forcingState->getField(snowName);
    rainAccessor.getHostRW() = forcingState->getField(rainName);

    // Postprocess the data fields
    auto& tair = tairAccessor.getAutoRW();
    auto& wind = windAccessor.getAutoRW();
    const auto& uWind = uwindAccessor.getAutoRO();
    const auto& vWind = vwindAccessor.getAutoRO();
    auto& lw = lw_inAccessor.getAutoRW();
    auto& sw = sw_inAccessor.getAutoRW();
    auto& rain = rainAccessor.getAutoRW();
    auto& snow = snowAccessor.getAutoRW();

    overElementsAuto(OVER_ELEMENTS_LAMBDA(const ElementIndex i) {
        // Temperatures in Celsius
        tair[i] -= Water::Tf;

        // Calculate wind speed from components
        wind[i] = std::sqrt(uWind[i] * uWind[i] + vWind[i] * vWind[i]);

        // Change accumulated values to instantaneous ones
        lw[i] /= 3._ft * 3600._ft;
        sw[i] /= 3._ft * 3600._ft;
        rain[i] /= 3._ft * 3600._ft;
        snow[i] /= 3._ft * 3600._ft;
    });

    fluxImpl->update(tst);

    try {
        checkFields();
    } catch (const std::exception& e) {
        throw std::runtime_error("CARRA2Atmosphere::update:: " + std::string(e.what()));
    }
}

void CARRA2Atmosphere::setFilePath(const std::string& filePathIn) { filePath = filePathIn; }

void CARRA2Atmosphere::setData(const ModelState::DataMap& ms)
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
