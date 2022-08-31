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

template <>
const std::map<int, std::string> Configured<ConfiguredOcean>::keyMap = {
    { ConfiguredOcean::SST_KEY, "ConfiguredOcean.sst" },
    { ConfiguredOcean::SSS_KEY, "ConfiguredOcean.sss" },
    { ConfiguredOcean::MLD_KEY, "ConfiguredOcean.mld" },
};

void ConfiguredOcean::configure()
{
    sst = Configured<ConfiguredOcean>::getConfiguration(
        Configured<ConfiguredOcean>::keyMap.at(SST_KEY), sst0);
    sss = Configured<ConfiguredOcean>::getConfiguration(
        Configured<ConfiguredOcean>::keyMap.at(SSS_KEY), sss0);
    mld = Configured<ConfiguredOcean>::getConfiguration(
        Configured<ConfiguredOcean>::keyMap.at(MLD_KEY), mld0);
}

} /* namespace Nextsim */
