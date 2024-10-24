/*!
 * @file ConfiguredOcean.cpp
 *
 * @date 24 Sep 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ConfiguredOcean.hpp"

#include "include/Finalizer.hpp"
#include "include/IFreezingPoint.hpp"
#include "include/IIceOceanHeatFlux.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/NextsimModule.hpp"
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
static const std::string uKey = pfx + ".current_u";
static const std::string vKey = pfx + ".current_v";

static const std::map<int, std::string> keyMap = {
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
    Module::getHelpRecursive<IIceOceanHeatFlux>(map, getAll);
    Module::getHelpRecursive<IFreezingPoint>(map, getAll);
    return map;
}

void ConfiguredOcean::configure()
{
    Finalizer::registerUnique(Module::finalize<IIceOceanHeatFlux>);
    Finalizer::registerUnique(Module::finalize<IFreezingPoint>);

    sst0 = Configured<ConfiguredOcean>::getConfiguration(keyMap.at(SST_KEY), sst0);
    sss0 = Configured<ConfiguredOcean>::getConfiguration(keyMap.at(SSS_KEY), sss0);
    mld0 = Configured<ConfiguredOcean>::getConfiguration(keyMap.at(MLD_KEY), mld0);
    u0 = Configured<ConfiguredOcean>::getConfiguration(keyMap.at(CURRENTU_KEY), u0);
    v0 = Configured<ConfiguredOcean>::getConfiguration(keyMap.at(CURRENTV_KEY), v0);

    // set the external SS* arrays as part of configuration, as opposed to at construction as normal
    getStore().registerArray(Protected::EXT_SST, &sstExt, RO);
    getStore().registerArray(Protected::EXT_SSS, &sssExt, RO);

    slabOcean.configure();

    tryConfigure(Module::getImplementation<IIceOceanHeatFlux>());
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

    /* It's only the SSH gradient which has an effect, so being able to sett a constant SSH is
     * useless. */
    ssh = 0.;

    slabOcean.setData(ms);

    Module::getImplementation<IIceOceanHeatFlux>().setData(ms);
}

void ConfiguredOcean::updateBefore(const TimestepTime& tst)
{
    Module::getImplementation<IIceOceanHeatFlux>().update(tst);
}

void ConfiguredOcean::updateAfter(const TimestepTime& tst)
{
    slabOcean.update(tst);
    sst = ModelArrayRef<Protected::SLAB_SST, RO>(getStore()).data();
    sss = ModelArrayRef<Protected::SLAB_SSS, RO>(getStore()).data();
}
} /* namespace Nextsim */
