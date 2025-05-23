/*!
 * @file TOPAZOcean.cpp
 *
 * @date 23 May 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/TOPAZOcean.hpp"

#include "include/Finalizer.hpp"
#include "include/NextsimModule.hpp"
#include "include/ParaGridIO.hpp"
#include "include/constants.hpp"
#include "include/gridNames.hpp"

namespace Nextsim {

std::string TOPAZOcean::filePath;

bool TOPAZOcean::checkFieldsDefault = false;

std::string TOPAZOcean::pfx = "TOPAZOcean";
std::string TOPAZOcean::fileKey = pfx + ".file";
std::string TOPAZOcean::checkFieldsKey = pfx + ".check_fields";

static const std::map<int, std::string> keyMap = {
    { TOPAZOcean::FILEPATH_KEY, TOPAZOcean::fileKey },
    { TOPAZOcean::CHECKFIELDS_KEY, TOPAZOcean::checkFieldsKey },
};

TOPAZOcean::TOPAZOcean()
    : sstExt(ModelArray::Type::H, { -5, 50 })
    , sssExt(ModelArray::Type::H, { 0, 50 })
    , slabOcean(m_couplingArrays)
{
}

ConfigurationHelp::HelpMap& TOPAZOcean::getHelpRecursive(HelpMap& map, bool getAll)
{
    map[pfx] = { { fileKey, ConfigType::STRING, {}, "", "",
                     "Path to the processed NetCDF file providing the TOPAZ forcings." },
        { checkFieldsKey, ConfigType::BOOLEAN, { "true", "false" },
            ConfigurationHelp::toString(checkFieldsDefault), "",
            "Set to true to check if the main module variables fall within a reasonable physical "
            "range." } };

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

    boolCheckFields = Configured::getConfiguration(keyMap.at(CHECKFIELDS_KEY), checkFieldsDefault);
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
    overElements(
        std::bind(&TOPAZOcean::updateTf, this, std::placeholders::_1, std::placeholders::_2),
        TimestepTime());

    pIOHeatFlux->update(tst);
}

void TOPAZOcean::updateAfter(const TimestepTime& tst)
{
    mergeFluxes(tst);
    slabOcean.update(tst);
    sst = ModelArrayRef<Protected::SLAB_SST, RO>(getStore());
    sss = ModelArrayRef<Protected::SLAB_SSS, RO>(getStore());

    if (boolCheckFields || boolCheckAll())
        checkFields(tst);
}

void TOPAZOcean::setFilePath(const std::string& filePathIn) { filePath = filePathIn; }

void TOPAZOcean::setData(const ModelState::DataMap& ms)
{
    IOceanBoundary::setData(ms);

    sstExt.resize();
    sssExt.resize();

    fieldsToCheck.emplace_back("sstExt", &sstExt);
    fieldsToCheck.emplace_back("sssExt", &sssExt);

    slabOcean.setData(ms);
}

void TOPAZOcean::updateTf(size_t i, const TimestepTime& tst) { tf[i] = (*pFreezingPoint)(sss[i]); }

} /* namespace Nextsim */
