/*!
 * @file SlabOcean.cpp
 *
 * @date 27 Jan 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/SlabOcean.hpp"

namespace Nextsim {

const std::string SlabOcean::sstSlabName = "sst_slab";
const std::string SlabOcean::sssSlabName = "sss_slab";
const double SlabOcean::defaultRelaxationTime = 30 * 24 * 60 * 60; // 30 days in seconds

// Configuration strings
static const std::string className = "SlabOcean";
static const std::string timeTName = "timeT";
static const std::string timeSName = "timeS";

template <>
const std::map<int, std::string> Configured<SlabOcean>::keyMap = {
    { SlabOcean::TIMET_KEY, className + "." + timeTName },
    { SlabOcean::TIMES_KEY, className + "." + timeSName },
};
void SlabOcean::configure()
{
    timeT = Configured::getConfiguration(keyMap.at(TIMET_KEY), defaultRelaxationTime);
    timeS = Configured::getConfiguration(keyMap.at(TIMES_KEY), defaultRelaxationTime);
}

SlabOcean::HelpMap& SlabOcean::getHelpText(HelpMap& map, bool getAll)
{
    map[className] = {
        { keyMap.at(TIMET_KEY), ConfigType::NUMERIC, { "0", "∞" },
            std::to_string(defaultRelaxationTime), "s",
            "Relaxation time of the slab ocean to external temperature forcing." },
        { keyMap.at(TIMES_KEY), ConfigType::NUMERIC, { "0", "∞" },
            std::to_string(defaultRelaxationTime), "s",
            "Relaxation time of the slab ocean to external salinity forcing." },
    };
    return map;
};

void SlabOcean::setData(const ModelState::DataMap& ms)
{
    sstSlab.resize();
    sssSlab.resize();
    qdw.resize();
    fdw.resize();
}

ModelState SlabOcean::getState() const
{
    return { {
                 { sstSlabName, sstSlab },
                 { sssSlabName, sssSlab },
             },
        {} };
}
ModelState SlabOcean::getState(const OutputLevel&) const { return getState(); }

std::unordered_set<std::string> SlabOcean::hFields() const { return { sstSlabName, sssSlabName }; }

void SlabOcean::update(const TimestepTime&)
{
    // Slab SST update
    // Slab SSS update
}

} /* namespace Nextsim */
