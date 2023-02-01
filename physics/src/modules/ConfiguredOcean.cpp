/*!
 * @file ConfiguredOcean.cpp
 *
 * @date Aug 31, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ConfiguredOcean.hpp"

#include "include/IFreezingPoint.hpp"
#include "include/IIceOceanHeatFlux.hpp"
#include "include/Module.hpp"
#include "include/constants.hpp"

namespace Nextsim {

double ConfiguredOcean::sst0 = -1.5;
double ConfiguredOcean::sss0 = 32;
double ConfiguredOcean::mld0 = 10;
double ConfiguredOcean::u0 = 0;
double ConfiguredOcean::v0 = 0;

static const std::string pfx = "ConfiguredOcean";
static const std::string sstKey = pfx + ".sst";
static const std::string sssKey = pfx + ".sss";
static const std::string mldKey = pfx + ".mld";
static const std::string uKey = pfx + "current_u";
static const std::string vKey = pfx + "current_v";

template <>
const std::map<int, std::string> Configured<ConfiguredOcean>::keyMap = {
    { ConfiguredOcean::SST_KEY, sstKey },
    { ConfiguredOcean::SSS_KEY, sssKey },
    { ConfiguredOcean::MLD_KEY, mldKey },
    { ConfiguredOcean::CURRENTU_KEY, uKey },
    { ConfiguredOcean::CURRENTV_KEY, vKey },
};

ConfiguredOcean::ConfiguredOcean()
    : sstExt(ModelArray::Type::H)
    , sssExt(ModelArray::Type::H)
{
}

ConfigurationHelp::HelpMap& ConfiguredOcean::getHelpRecursive(HelpMap& map, bool getAll)
{
    map[pfx] = {
        { sstKey, ConfigType::NUMERIC, { "-273", "374" }, std::to_string(sst0), "",
            "Sea surface temperature (˚C)." },
        { sssKey, ConfigType::NUMERIC, { "0", "1000" }, std::to_string(sss0), "",
            "Sea surface salinity (PSU)." },
        { mldKey, ConfigType::NUMERIC, { "0", "10984" }, std::to_string(mld0), "",
            "Mixed layer depth (m)." },
        { uKey, ConfigType::NUMERIC, { "-∞", "∞" }, std::to_string(u0), "",
            "Ocean current in the x (eastward) direction (m s⁻¹)." },
        { vKey, ConfigType::NUMERIC, { "-∞", "∞" }, std::to_string(v0), "",
            "Ocean current in the y (northward) direction (m s⁻¹)." },

    };
    return map;
}

void ConfiguredOcean::configure()
{
    sst0 = Configured<ConfiguredOcean>::getConfiguration(
        Configured<ConfiguredOcean>::keyMap.at(SST_KEY), sst0);
    sss0 = Configured<ConfiguredOcean>::getConfiguration(
        Configured<ConfiguredOcean>::keyMap.at(SSS_KEY), sss0);
    mld0 = Configured<ConfiguredOcean>::getConfiguration(
        Configured<ConfiguredOcean>::keyMap.at(MLD_KEY), mld0);
    u0 = Configured<ConfiguredOcean>::getConfiguration(
        Configured<ConfiguredOcean>::keyMap.at(CURRENTU_KEY), u0);
    v0 = Configured<ConfiguredOcean>::getConfiguration(
        Configured<ConfiguredOcean>::keyMap.at(CURRENTV_KEY), v0);

    // set the external SS* arrays as part of configuration, as opposed to at construction as normal
    registerProtectedArray(ProtectedArray::EXT_SST, &sstExt);
    registerProtectedArray(ProtectedArray::EXT_SSS, &sssExt);

    slabOcean.configure();
}

void ConfiguredOcean::setData(const ModelState::DataMap& ms)
{
    IOceanBoundary::setData(ms);
    sstExt.resize();
    sssExt.resize();

    sstExt = sst0;
    sssExt = sss0;
    mld = mld0;
    u = u0;
    v = v0;
    tf = Module::getImplementation<IFreezingPoint>()(sssExt[0]);
    cpml = Water::rho * Water::cp * mld[0];

    slabOcean.setData(ms);
}

void ConfiguredOcean::updateBefore(const TimestepTime& tst)
{
    Module::getImplementation<IIceOceanHeatFlux>().update(tst);
}

void ConfiguredOcean::updateAfter(const TimestepTime& tst)
{
    slabOcean.update(tst);
    sst = *getProtectedArray()[static_cast<size_t>(ProtectedArray::SLAB_SST)];
    sss = *getProtectedArray()[static_cast<size_t>(ProtectedArray::SLAB_SSS)];
}
} /* namespace Nextsim */
