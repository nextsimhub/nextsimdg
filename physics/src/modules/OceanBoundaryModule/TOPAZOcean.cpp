/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "include/TOPAZOcean.hpp"

#include "include/Finalizer.hpp"
#include "include/NextsimModule.hpp"
#include "include/ParaGridIO.hpp"
#include "include/constants.hpp"
#include "include/gridNames.hpp"

// testing
#include "kokkos/include/KokkosTimer.hpp"

namespace Nextsim {

std::string TOPAZOcean::filePath;

static const std::string pfx = "TOPAZOcean";
static const std::string fileKey = pfx + ".file";

static const std::map<int, std::string> keyMap = {
    { TOPAZOcean::FILEPATH_KEY, fileKey },
};

TOPAZOcean::TOPAZOcean()
    : sstExtAccessor(getStore(), RO, ModelArray::Type::H, std::pair(-5.0, 50.0))
    , sssExtAccessor(getStore(), RO, ModelArray::Type::H, std::pair(0.0, 50.0))
    , sstSlabAccessor(getStore())
    , sssSlabAccessor(getStore())
    , slabOcean(m_couplingArrays)
{
}

ConfigurationHelp::HelpMap& TOPAZOcean::getHelpRecursive(HelpMap& map, bool getAll)
{
    map[pfx] = {
        { fileKey, ConfigType::STRING, {}, "", "",
            "Path to the processed NetCDF file providing the TOPAZ forcings." },
    };

    return map;
}

void TOPAZOcean::configure()
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

ConfigMap TOPAZOcean::getConfiguration() const { return { { keyMap.at(FILEPATH_KEY), filePath } }; }

static const std::set<std::string> forcings = { sstName, sssName, mldName, uName, vName, sshName };

void TOPAZOcean::updateBefore(const TimestepTime& tst)
{
    // Read TOPAZ forcings at midnight
    if (std::fmod((tst.start - TimePoint()).seconds(), 86400.) == 0.) {
        forcingState = ParaGridIO::readForcingTimeStatic(forcings, tst.start, filePath);
    }

    sstExtAccessor.getHostRW() = forcingState.data.at(sstName);
    sssExtAccessor.getHostRW() = forcingState.data.at(sssName);
    mldAccessor.getHostRW() = forcingState.data.at(mldName);
    uAccessor.getHostRW() = forcingState.data.at(uName);
    vAccessor.getHostRW() = forcingState.data.at(vName);
    HField& ssh = sshAccessor.getHostRW();

    if (forcingState.data.count(sshName)) {
        ssh = forcingState.data.at(sshName);
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

void TOPAZOcean::updateAfter(const TimestepTime& tst)
{
    mergeFluxes(tst);
    slabOcean.update(tst);

    sstAccessor.getAutoRW().assignData(sstSlabAccessor.getAutoRO());
    sssAccessor.getAutoRW().assignData(sssSlabAccessor.getAutoRO());

    try {
        checkFields();
    } catch (const std::exception& e) {
        throw std::runtime_error("TOPAZOcean::updateAfter: " + std::string(e.what()));
    }
}

ModelState TOPAZOcean::getStatePrognostic() const
{
    ModelState state = IOceanBoundary::getStatePrognostic();
    return state.merge(slabOcean.getStatePrognostic());
}

ModelState TOPAZOcean::getStateDiagnostic() const
{
    ModelState state = IOceanBoundary::getStateDiagnostic();
    return state.merge(slabOcean.getStateDiagnostic());
}

void TOPAZOcean::setFilePath(const std::string& filePathIn) { filePath = filePathIn; }

void TOPAZOcean::setData(const ModelState::DataMap& ms)
{
    IOceanBoundary::setData(ms);

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

// void TOPAZOcean::updateTf(size_t i, const TimestepTime& tst) { tf[i] = (*pFreezingPoint)(sss[i]);
// }

} /* namespace Nextsim */
