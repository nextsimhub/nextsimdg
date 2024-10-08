/*!
 * @file TOPAZOcean.cpp
 *
 * @date 24 Sep 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/TOPAZOcean.hpp"

#include "include/IFreezingPoint.hpp"
#include "include/IIceOceanHeatFlux.hpp"
#include "include/Module.hpp"
#include "include/ParaGridIO.hpp"
#include "include/constants.hpp"

namespace Nextsim {

std::string TOPAZOcean::filePath;

static const std::string pfx = "TOPAZOcean";
static const std::string fileKey = pfx + ".file";

template <>
const std::map<int, std::string> Configured<TOPAZOcean>::keyMap = {
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
    filePath = Configured::getConfiguration(keyMap.at(FILEPATH_KEY), std::string());

    slabOcean.configure();

    getStore().registerArray(Protected::EXT_SST, &sstExt, RO);
    getStore().registerArray(Protected::EXT_SSS, &sssExt, RO);
}

void TOPAZOcean::updateBefore(const TimestepTime& tst)
{
    // TODO: Get more authoritative names for the forcings
    std::set<std::string> forcings = { "sst", "sss", "mld", "u", "v", "ssh" };

    ModelState state = ParaGridIO::readForcingTimeStatic(forcings, tst.start, filePath);
    sstExt = state.data.at("sst");
    sssExt = state.data.at("sss");
    mld = state.data.at("mld");
    u = state.data.at("u");
    v = state.data.at("v");
    ssh = state.data.at("ssh");

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
