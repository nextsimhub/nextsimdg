/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#include "include/TOPAZOcean.hpp"

#include "include/Finalizer.hpp"
#include "include/NextsimModule.hpp"
#include "include/ParaGridIO.hpp"
#include "include/constants.hpp"
#include "include/gridNames.hpp"

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

void TOPAZOcean::updateBefore(const TimestepTime& tst)
{
    std::set<std::string> forcings = { sstName, sssName, mldName, uName, vName, sshName };

    ModelState state = ParaGridIO::readForcingTimeStatic(forcings, tst.start, filePath);
    sstExtAccessor.getHostRW() = state.data.at(sstName);
    sssExtAccessor.getHostRW() = state.data.at(sssName);
    HField& mld = mldAccessor.getHostRW();
    mld = state.data.at(mldName);
    uAccessor.getHostRW() = state.data.at(uName);
    vAccessor.getHostRW() = state.data.at(vName);
    HField& ssh = sshAccessor.getHostRW();
    if (state.data.count(sshName)) {
        ssh = state.data.at(sshName);
    } else {
        ssh = 0.;
    }

    cpmlAccessor.getHostRW() = Water::rhoOcean * Water::cp * mld;

    // Update the freezing point
    const HField& sss = sssAccessor.getHostRO();
    HField& tf = tfAccessor.getHostRW();
    overElements(
        [&](size_t i, const TimestepTime& tsTime) { tf[i] = (*pFreezingPoint)(sss[i]); }, tst);

    Module::getImplementation<IIceOceanHeatFlux>().update(tst);
}

void TOPAZOcean::updateAfter(const TimestepTime& tst)
{
    mergeFluxes(tst);
    slabOcean.update(tst);
    sstAccessor.getHostRW().component() = sstSlabAccessor.getHostRO().component();
    sssAccessor.getHostRW().component() = sssSlabAccessor.getHostRO().component();

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
    sstExt.resize();
    HField& sssExt = sssExtAccessor.getHostRW();
    sssExt.resize();

    addChecks({
        { "sstExt", &sstExt },
        { "sssExt", &sssExt },
    });

    slabOcean.setData(ms);
}

// void TOPAZOcean::updateTf(size_t i, const TimestepTime& tst) { tf[i] = (*pFreezingPoint)(sss[i]);
// }

} /* namespace Nextsim */
