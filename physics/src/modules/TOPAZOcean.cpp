/*!
 * @file TOPAZOcean.cpp
 *
 * @date Nov 25, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/TOPAZOcean.hpp"

#include "include/Module.hpp"
#include "include/ParaGridIO.hpp"

#include <iostream> // FIXME remove me
namespace Nextsim {

std::string TOPAZOcean::filePath;

static const std::string pfx = "TOPAZOcean";
static const std::string fileKey = pfx + ".file";

template <>
const std::map<int, std::string> Configured<TOPAZOcean>::keyMap = {
    { TOPAZOcean::FILEPATH_KEY, fileKey },
};

TOPAZOcean::TOPAZOcean()
{
    registerProtectedArray(ProtectedArray::SST, &sst);
    registerProtectedArray(ProtectedArray::SSS, &sss);
    registerProtectedArray(ProtectedArray::MLD, &mld);
    registerProtectedArray(ProtectedArray::OCEAN_U, &u);
    registerProtectedArray(ProtectedArray::OCEAN_V, &v);
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
}

void TOPAZOcean::updateBefore(const TimestepTime& tst)
{
    // TODO: Get more authoritative names for the forcings
    std::set<std::string> forcings = { "sst", "sss", "mld", "u", "v" };

    ModelState state = ParaGridIO::readForcingTimeStatic(forcings, tst.start, filePath);
    sst = state.data.at("sst");
    sss = state.data.at("sss");
    mld = state.data.at("mld");
    u = state.data.at("u");
    v = state.data.at("v");
}

void TOPAZOcean::setFilePath(const std::string& filePathIn) { filePath = filePathIn; }

void TOPAZOcean::setData(const ModelState::DataMap&) { }

} /* namespace Nextsim */
