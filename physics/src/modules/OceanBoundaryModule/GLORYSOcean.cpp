/*!
 * @author  Einar Ólason <einar.olason@nersc.no>
 */

#include "include/GLORYSOcean.hpp"

#include "include/Finalizer.hpp"
#include "include/NextsimModule.hpp"
#include "include/ParaGridIO.hpp"
#include "include/constants.hpp"

namespace Nextsim {

std::string GLORYSOcean::filePath;

static const std::string pfx = "GLORYSOcean";
static const std::string fileKey = pfx + ".file";

static const std::map<int, std::string> keyMap = {
    { GLORYSOcean::FILEPATH_KEY, fileKey },
};

GLORYSOcean::GLORYSOcean()
    : sstExtAccessor(getStore(), RO, ModelArray::Type::H, std::pair(-5.0, 50.0))
    , sssExtAccessor(getStore(), RO, ModelArray::Type::H, std::pair(0.0, 50.0))
    , sstSlabAccessor(getStore())
    , sssSlabAccessor(getStore())
    , slabOcean(m_couplingArrays)
{
}

ConfigurationHelp::HelpMap& GLORYSOcean::getHelpRecursive(HelpMap& map, bool getAll)
{
    map[pfx] = {
        { fileKey, ConfigType::STRING, {}, "", "",
            "Path to the processed NetCDF file providing the GLORYS forcings." },
    };

    return map;
}

void GLORYSOcean::configure()
{
    Finalizer::registerUnique(Module::finalize<IIceOceanHeatFlux>);
    Finalizer::registerUnique(Module::finalize<IFreezingPoint>);

    pIOHeatFlux = std::move(Module::getInstance<IIceOceanHeatFlux>());
    tryConfigure(*pIOHeatFlux);

    pFreezingPoint = std::move(Module::getInstance<IFreezingPoint>());
    tryConfigure(*pFreezingPoint);

    filePath = Configured::getConfiguration(keyMap.at(FILEPATH_KEY), std::string());

    slabOcean.configure();
}

ConfigMap GLORYSOcean::getConfiguration() const
{
    return { { keyMap.at(FILEPATH_KEY), filePath } };
}

void GLORYSOcean::updateBefore(const TimestepTime& tst)
{
    forcingState->update(tst.start);

    sstExtAccessor.getHostRW() = forcingState->getField(sstName);
    sssExtAccessor.getHostRW() = forcingState->getField(sssName);
    mldAccessor.getHostRW() = forcingState->getField(mldName);
    uAccessor.getHostRW() = forcingState->getField(uName);
    vAccessor.getHostRW() = forcingState->getField(vName);
    HField& ssh = sshAccessor.getHostRW();

    if (forcingState->isValid(sshName)) {
        ssh = forcingState->getField(sshName);
    } else {
        ssh = 0.;
    }

    auto& cpml = cpmlAccessor.getAutoRW();
    const auto& mld = mldAccessor.getAutoRO();
    overElementsAuto(OVER_ELEMENTS_LAMBDA(
        const ElementIndex i) { cpml[i] = Water::rhoOcean * Water::cp * mld[i]; });

    // Update the freezing point
    auto& tf = tfAccessor.getAutoRW();
    const auto& sss = sssAccessor.getAutoRO();
    pFreezingPoint->update(tf, sss);

    Module::getImplementation<IIceOceanHeatFlux>().update(tst);
}

void GLORYSOcean::updateAfter(const TimestepTime& tst)
{
    mergeFluxes(tst);
    slabOcean.update(tst);

    sstAccessor.getAutoRW().assignData(sstSlabAccessor.getAutoRO());
    sssAccessor.getAutoRW().assignData(sssSlabAccessor.getAutoRO());

    try {
        checkFields();
    } catch (const std::exception& e) {
        throw std::runtime_error("GLORYSOcean::updateAfter: " + std::string(e.what()));
    }
}

ModelState GLORYSOcean::getStatePrognostic() const
{
    ModelState state = IOceanBoundary::getStatePrognostic();
    return state.merge(slabOcean.getStatePrognostic());
}

ModelState GLORYSOcean::getStateDiagnostic() const
{
    ModelState state = IOceanBoundary::getStateDiagnostic();
    return state.merge(slabOcean.getStateDiagnostic());
}

void GLORYSOcean::setFilePath(const std::string& filePathIn) { filePath = filePathIn; }

void GLORYSOcean::setData(const ModelState::DataMap& ms)
{
    IOceanBoundary::setData(ms);

    ModelState state;
    const ModelMetadata& metadata = ModelMetadata::getInstance();
    metadata.affixCoordinates(state);
    forcingState->initialise(metadata.startTime(), filePath, ncLonName, ncLatName, ncTimeName,
        forcings, vectors, state.data.at(longitudeName), state.data.at(latitudeName));

    HField& sstExt = sstExtAccessor.getHostRW();
    sstExt.reinitialize();
    HField& sssExt = sssExtAccessor.getHostRW();
    sssExt.reinitialize();

    addChecks({
        { "sstExt", sstExtAccessor },
        { "sssExt", sssExtAccessor },
    });

    slabOcean.setData(ms);
}

} /* namespace Nextsim */
