/*!
 * @file TOPAZOcean.cpp
 *
 * @date 17 Oct 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/TOPAZOcean.hpp"

#include "include/Finalizer.hpp"
#include "include/IIceOceanHeatFlux.hpp"
#include "include/IFreezingPoint.hpp"
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
    : sstExt(ModelArray::Type::H)
    , sssExt(ModelArray::Type::H)
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

    filePath = Configured::getConfiguration(keyMap.at(FILEPATH_KEY), std::string());

    slabOcean.configure();

    getStore().registerArray(Protected::EXT_SST, &sstExt, RO);
    getStore().registerArray(Protected::EXT_SSS, &sssExt, RO);
}

void TOPAZOcean::updateBefore(const TimestepTime& tst)
{
    // TODO: Get more authoritative names for the forcings
    std::set<std::string> forcings = { sstName, sssName, "mld", uName, vName, sshName };

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

    cpml = Water::rho * Water::cp * mld;
    overElements(
        std::bind(&TOPAZOcean::updateTf, this, std::placeholders::_1, std::placeholders::_2),
        TimestepTime());

    Module::getImplementation<IIceOceanHeatFlux>().update(tst);
}

void TOPAZOcean::updateAfter(const TimestepTime& tst)
{
    slabOcean.update(tst);
    sst = ModelArrayRef<Protected::SLAB_SST, RO>(getStore()).data();
    sss = ModelArrayRef<Protected::SLAB_SSS, RO>(getStore()).data();
}

void TOPAZOcean::setFilePath(const std::string& filePathIn) { filePath = filePathIn; }

void TOPAZOcean::setData(const ModelState::DataMap& ms)
{
    IOceanBoundary::setData(ms);

    sstExt.resize();
    sssExt.resize();
    slabOcean.setData(ms);
}

void TOPAZOcean::updateTf(size_t i, const TimestepTime& tst)
{
    tf[i] = Module::getImplementation<IFreezingPoint>()(sss[i]);
}

} /* namespace Nextsim */
