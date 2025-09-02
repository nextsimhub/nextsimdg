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
    : sstExt(ModelArray::Type::H, { -5, 50 })
    , sssExt(ModelArray::Type::H, { 0, 50 })
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

    getStore().registerArray(Protected::EXT_SST, &sstExt, RO);
    getStore().registerArray(Protected::EXT_SSS, &sssExt, RO);
}

ConfigMap TOPAZOcean::getConfiguration() const { return { { keyMap.at(FILEPATH_KEY), filePath } }; }

void TOPAZOcean::updateBefore(const TimestepTime& tst)
{
    std::set<std::string> forcings = { sstName, sssName, mldName, uName, vName, sshName };

    ModelState state = ParaGridIO::readForcingTimeStatic(forcings, tst.start, filePath);
    sstExt = state.data.at(sstName);
    sssExt = state.data.at(sssName);
    mld = state.data.at(mldName);
    u = state.data.at(uName);
    v = state.data.at(vName);
    if (state.data.count(sshName)) {
        ssh = state.data.at(sshName);
    } else {
        ssh = 0.;
    }

    cpml = Water::rhoOcean * Water::cp * mld;
    overElements([this](size_t i, const TimestepTime& tsTime) { this->updateTf(i, tsTime); }, tst);

    Module::getImplementation<IIceOceanHeatFlux>().update(tst);
}

void TOPAZOcean::updateAfter(const TimestepTime& tst)
{
    mergeFluxes(tst);
    slabOcean.update(tst);
    sst = ModelArrayRef<Protected::SLAB_SST, RO>(getStore());
    sss = ModelArrayRef<Protected::SLAB_SSS, RO>(getStore());

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

    sstExt.resize();
    sssExt.resize();

    addChecks({
        { "sstExt", &sstExt },
        { "sssExt", &sssExt },
    });

    slabOcean.setData(ms);

    // We assume that all incoming shortwave radiation is absorbed in the mixed layer
    fracQSWAbs = 1.;
}

void TOPAZOcean::updateTf(size_t i, const TimestepTime& tst) { tf[i] = (*pFreezingPoint)(sss[i]); }

} /* namespace Nextsim */
