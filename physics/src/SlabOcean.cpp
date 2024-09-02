/*!
 * @file SlabOcean.cpp
 *
 * @date 30 Aug 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/SlabOcean.hpp"

#include "include/constants.hpp"

#include <map>
#include <string>

namespace Nextsim {

const std::string SlabOcean::sstSlabName = "sst_slab";
const std::string SlabOcean::sssSlabName = "sss_slab";
const double SlabOcean::defaultRelaxationTime = 30 * 24 * 60 * 60; // 30 days in seconds

// Configuration strings
static const std::string className = "SlabOcean";
static const std::string relaxationTimeTName = "timeT";
static const std::string relaxationTimeSName = "timeS";

template <>
const std::map<int, std::string> Configured<SlabOcean>::keyMap = {
    { SlabOcean::TIMET_KEY, className + "." + relaxationTimeTName },
    { SlabOcean::TIMES_KEY, className + "." + relaxationTimeSName },
};
void SlabOcean::configure()
{
    relaxationTimeT = Configured::getConfiguration(keyMap.at(TIMET_KEY), defaultRelaxationTime);
    relaxationTimeS = Configured::getConfiguration(keyMap.at(TIMES_KEY), defaultRelaxationTime);

    getStore().registerArray(Protected::SLAB_QDW, &qdw, RO);
    getStore().registerArray(Protected::SLAB_FDW, &fdw, RO);
    getStore().registerArray(Protected::SLAB_SST, &sstSlab, RO);
    getStore().registerArray(Protected::SLAB_SSS, &sssSlab, RO);
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
    qdw.resize();
    fdw.resize();
    sstSlab.resize();
    sssSlab.resize();
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

void SlabOcean::update(const TimestepTime& tst)
{
    double dt = tst.step.seconds();

    // Slab SST update
    qdw = (sstExt - sst) * cpml / relaxationTimeT;
    HField qioMean = qio * cice; // cice at start of TS, not updated
    HField qowMean = qow * (1 - cice); // 1- cice = open water fraction
    sstSlab = sst - dt * (qioMean + qowMean - qdw) / cpml;

    // Slab SSS update
    HField arealDensity = cpml / Water::cp; // density times depth, or cpml divided by cp
    // This is simplified compared to the finiteelement.cpp calculation
    // Fdw = delS * mld * physical::rhow /(timeS*M_sss[i] - ddt*delS) where delS = sssSlab - sssExt
    fdw = (1 - sssExt / sss) * arealDensity / relaxationTimeS;

    /* We use the "virtual salinity flux" approach:
     *   \partial S / \partial t = F_{fw} S / H,
     * where S is the ocean salinity, F_{fw} is the freshwater flux and H the mixed-layer depth. In
     * the presence of sea ice F_{fw} S must be replaced by (S-S_i)I + S(P-E+M_s), with I and S
     * the fresh-water flux due to ice and snow melt, respectively, and P and E, as precipitation
     * and evaporation, respectively. Using the salt (F_s) and fresh-water fluxes (F_{fw})
     * calculated by the IOceanBoundary class, we can write the change in ocean salinity as
     *   \partial S / \partial t = ( S F_{fw} - 10^3 F_s ) / H
     */
    sssSlab = sss + (sss * (fwFlux + fdw) - 1e3 * sFlux) * dt / arealDensity;
}

} /* namespace Nextsim */
