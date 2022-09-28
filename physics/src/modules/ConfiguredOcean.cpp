/*!
 * @file ConfiguredOcean.cpp
 *
 * @date Aug 31, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ConfiguredOcean.hpp"

namespace Nextsim {

double ConfiguredOcean::sst0 = -1.5;
double ConfiguredOcean::sss0 = 32;
double ConfiguredOcean::mld0 = 10;
double ConfiguredOcean::u0 = 0;
double ConfiguredOcean::v0 = 0;

template <>
const std::map<int, std::string> Configured<ConfiguredOcean>::keyMap = {
    { ConfiguredOcean::SST_KEY, "ConfiguredOcean.sst" },
    { ConfiguredOcean::SSS_KEY, "ConfiguredOcean.sss" },
    { ConfiguredOcean::MLD_KEY, "ConfiguredOcean.mld" },
    { ConfiguredOcean::CURRENTU_KEY, "ConfiguredOcean.current_u" },
    { ConfiguredOcean::CURRENTV_KEY, "ConfiguredOcean.current_v" },
};

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
}

void ConfiguredOcean::setData(const ModelState::DataMap& ms)
{
    IOceanBoundary::setData(ms);
    sst = sst0;
    sss = sss0;
    mld = mld0;
    u = u0;
    v = v0;
}

} /* namespace Nextsim */
