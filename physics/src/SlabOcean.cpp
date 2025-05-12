/*!
 * @file SlabOcean.cpp
 *
 * @date 29 Apr 2025
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

static const std::map<int, std::string> keyMap = {
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
            ConfigurationHelp::toString(defaultRelaxationTime), "s",
            "Relaxation time of the slab ocean to external temperature forcing." },
        { keyMap.at(TIMES_KEY), ConfigType::NUMERIC, { "0", "∞" },
            ConfigurationHelp::toString(defaultRelaxationTime), "s",
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

void SlabOcean::update(const TimestepTime& tst)
{
    dt = tst.step.seconds();
    overElements(
        std::bind(&SlabOcean::updateElement, this, std::placeholders::_1, std::placeholders::_2),
        tst);
}

void SlabOcean::updateElement(size_t i, const TimestepTime& tst)
{
    // Slab SST update
    qdw[i] = (sstExt[i] - sst[i]) * cpml[i] / relaxationTimeT;
    sstSlab[i] = sst[i] - dt * (qswNet[i] + qNoSun[i] - qdw[i]) / cpml[i];

    // Slab SSS update
    const double arealDensity = cpml[i] / Water::cp; // density times depth, or cpml divided by cp
    // This is simplified compared to the finiteelement.cpp calculation
    // Fdw = delS * mld * physical::rhow /(timeS*M_sss[i] - ddt*delS) where delS = sssSlab - sssExt
    fdw[i] = (1 - sssExt[i] / sss[i]) * arealDensity / relaxationTimeS;

    // Mass per unit area after all the changes in water volume
    // Clamp the denominator to be at least 1 m deep, i.e. at least Water::rho kg m⁻²
    const double denominator = std::max(arealDensity - (fwFlux[i] - fdw[i]) * dt, Water::rhoOcean);
    sssSlab[i] = sss[i] + (sss[i] * fwFlux[i] - fdw[i] * dt) / denominator;
}

} /* namespace Nextsim */
