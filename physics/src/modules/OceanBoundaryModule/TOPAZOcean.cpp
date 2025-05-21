/*!
 * @file TOPAZOcean.cpp
 *
 * @date 12 Mar 2025
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
std::string TOPAZOcean::fieldNamesKey = pfx + ".fields_names";
std::string TOPAZOcean::fieldNamesDefault = "sst_ext,sss_ext,mld,ocean_u,ocean_v,ssh";

static const std::map<int, std::string> keyMap = {
    { TOPAZOcean::FILEPATH_KEY, TOPAZOcean::fileKey },
    { TOPAZOcean::FIELDNAMES_KEY, TOPAZOcean::fieldNamesKey },
    { TOPAZOcean::CHECKFIELDS_KEY, TOPAZOcean::checkFieldsKey },
};

TOPAZOcean::TOPAZOcean()
    : sstExt(ModelArray::Type::H)
    , sssExt(ModelArray::Type::H)
    , slabOcean(m_couplingArrays)
{
}

ConfigurationHelp::HelpMap& TOPAZOcean::getHelpRecursive(HelpMap& map, bool getAll)
{
    auto options
        = getCheckingHelpList(fieldNamesKey, fieldNamesDefault, checkFieldsKey, checkFieldsDefault);
    options.push_back({ fileKey, ConfigType::STRING, {}, "", "",
        "Path to the processed NetCDF file providing the TOPAZ forcings." });
    map[pfx] = options;

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

    if (getConfiguration(keyMap.at(CHECKFIELDS_KEY), checkFieldsDefault)) {
        setFieldsToCheck(fieldNamesDefault, pfx);
    }
}

ConfigMap TOPAZOcean::getConfiguration() const
{
    return { { keyMap.at(FILEPATH_KEY), filePath } };
}

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

    checkFields(tst);
}

void TOPAZOcean::updateAfter(const TimestepTime& tst)
{
    mergeFluxes(tst);
    slabOcean.update(tst);
    sst = ModelArrayRef<Protected::SLAB_SST, RO>(getStore());
    sss = ModelArrayRef<Protected::SLAB_SSS, RO>(getStore());
}

void TOPAZOcean::setFilePath(const std::string& filePathIn) { filePath = filePathIn; }

void TOPAZOcean::setData(const ModelState::DataMap& ms)
{
    IOceanBoundary::setData(ms);

    sstExt.resize();
    sssExt.resize();
    slabOcean.setData(ms);
}

void TOPAZOcean::updateTf(size_t i, const TimestepTime& tst) { tf[i] = (*pFreezingPoint)(sss[i]); }

} /* namespace Nextsim */
