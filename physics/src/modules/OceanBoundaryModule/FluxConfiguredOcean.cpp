/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#include "include/FluxConfiguredOcean.hpp"

#include "include/IFreezingPoint.hpp"
#include "include/NextsimModule.hpp"
#include "include/constants.hpp"

namespace Nextsim {

FloatType FluxConfiguredOcean::qio0 = 0;
FloatType FluxConfiguredOcean::sst0 = -1.5;
FloatType FluxConfiguredOcean::sss0 = 32;
FloatType FluxConfiguredOcean::mld0 = 10;
FloatType FluxConfiguredOcean::u0 = 0;
FloatType FluxConfiguredOcean::v0 = 0;

static const std::string pfx = "FluxConfiguredOcean";
static const std::string qioKey = pfx + ".qio";
static const std::string sstKey = pfx + ".sst";
static const std::string sssKey = pfx + ".sss";
static const std::string mldKey = pfx + ".mld";
static const std::string uKey = pfx + ".current_u";
static const std::string vKey = pfx + ".current_v";

static const std::map<int, std::string> keyMap = {
    { FluxConfiguredOcean::QIO_KEY, qioKey },
    { FluxConfiguredOcean::SST_KEY, sstKey },
    { FluxConfiguredOcean::SSS_KEY, sssKey },
    { FluxConfiguredOcean::MLD_KEY, mldKey },
    { FluxConfiguredOcean::CURRENTU_KEY, uKey },
    { FluxConfiguredOcean::CURRENTV_KEY, vKey },
};

ConfigurationHelp::HelpMap& FluxConfiguredOcean::getHelpRecursive(HelpMap& map, bool getAll)
{
    map[pfx] = {
        { qioKey, ConfigType::NUMERIC, { "-∞", "∞" }, ConfigurationHelp::toString(qio0), "",
            "Ocean to ice heat flux (W m⁻²)." },
        { sstKey, ConfigType::NUMERIC, { "-273", "374" }, ConfigurationHelp::toString(sst0), "",
            "Sea surface temperature (˚C)." },
        { sssKey, ConfigType::NUMERIC, { "0", "1000" }, ConfigurationHelp::toString(sss0), "",
            "Sea surface salinity (PSU)." },
        { mldKey, ConfigType::NUMERIC, { "0", "10984" }, ConfigurationHelp::toString(mld0), "",
            "Mixed layer depth (m)." },
        { uKey, ConfigType::NUMERIC, { "-∞", "∞" }, ConfigurationHelp::toString(u0), "",
            "Ocean current in the x (eastward) direction (m s⁻¹)." },
        { vKey, ConfigType::NUMERIC, { "-∞", "∞" }, ConfigurationHelp::toString(v0), "",
            "Ocean current in the y (northward) direction (m s⁻¹)." },

    };
    return map;
}

void FluxConfiguredOcean::configure()
{
    qio0 = Configured<FluxConfiguredOcean>::getConfiguration(keyMap.at(QIO_KEY), qio0);
    sst0 = Configured<FluxConfiguredOcean>::getConfiguration(keyMap.at(SST_KEY), sst0);
    sss0 = Configured<FluxConfiguredOcean>::getConfiguration(keyMap.at(SSS_KEY), sss0);
    mld0 = Configured<FluxConfiguredOcean>::getConfiguration(keyMap.at(MLD_KEY), mld0);
    u0 = Configured<FluxConfiguredOcean>::getConfiguration(keyMap.at(CURRENTU_KEY), u0);
    v0 = Configured<FluxConfiguredOcean>::getConfiguration(keyMap.at(CURRENTV_KEY), v0);
}

ConfigMap FluxConfiguredOcean::getConfiguration() const
{
    return {
        { keyMap.at(QIO_KEY), qio0 },
        { keyMap.at(SST_KEY), sst0 },
        { keyMap.at(SSS_KEY), sss0 },
        { keyMap.at(MLD_KEY), mld0 },
        { keyMap.at(CURRENTU_KEY), u0 },
        { keyMap.at(CURRENTV_KEY), v0 },
    };
}
void FluxConfiguredOcean::setData(const ModelState::DataMap& ms)
{
    IOceanBoundary::setData(ms);
    qioAccessor.getHostRW() = qio0;
    sstAccessor.getHostRW() = sst0;
    HField& sss = sssAccessor.getHostRW();
    sss = sss0;
    HField& mld = mldAccessor.getHostRW();
    mld = mld0;
    uAccessor.getHostRW() = u0;
    vAccessor.getHostRW() = v0;
    tfAccessor.getHostRW() = Module::getImplementation<IFreezingPoint>()(sss[0]);
    cpmlAccessor.getHostRW() = Water::rho * Water::cp * mld[0];

    /* It's only the SSH gradient which has an effect, so being able to set a constant SSH is
     * useless. */
    sshAccessor.getHostRW() = 0.;
}

} /* namespace Nextsim */
