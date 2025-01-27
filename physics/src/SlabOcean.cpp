/*!
 * @file SlabOcean.cpp
 *
 * @date 27 Jan 2025
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

    getStore().registerArray(Protected::SLAB_FDW, &fdw, RO);
    getStore().registerArray(Protected::SLAB_QDW, &qdw, RO);
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
}

ModelState SlabOcean::getState() const { return { {}, { { getConfiguration() } } }; }

ModelState SlabOcean::getState(const OutputLevel&) const { return getState(); }

void SlabOcean::update(const TimestepTime& tst)
{
    overElements(
        std::bind(&SlabOcean::updateElements, this, std::placeholders::_1, std::placeholders::_2),
        TimestepTime());
}

void SlabOcean::updateElements(size_t i, const TimestepTime& tst)
{
    double dt = tst.step.seconds();

    // Slab SST update
    qdw[i] = (sstExt[i] - sst[i]) * cpml[i] / relaxationTimeT;
    const double qioMean = qio[i] * cice[i]; // cice at start of TS, not updated
    const double qowMean = qow[i] * (1 - cice[i]); // 1- cice = open water fraction
    sst[i] -= dt * (qioMean + qowMean - qdw[i]) / cpml[i];

    // Slab SSS update
    const double arealDensity = cpml[i] / Water::cp; // density times depth, or cpml divided by cp

    // This is simplified compared to the finiteelement.cpp calculation
    // Fdw = delS * mld * physical::rhow /(timeS*M_sss[i] - ddt*delS) where delS = sssSlab - sssExt
    fdw[i] = (1 - sssExt[i] / sss[i]) * arealDensity / relaxationTimeS;

    // ice volume change, both laterally and vertically
    const double deltaIceVol = newIce[i] + deltaHice[i] * cice[i];

    // change in snow volume due to melting (should be < 0)
    const double meltSnowVol = deltaSmelt[i] * cice[i];

    // Mass per unit area after all the changes in water volume
    // Clamp the denominator to be at least 1 m deep, i.e. at least Water::rho kg m⁻²
    const double denominator = std::max(Water::rho,
        arealDensity - deltaIceVol * Ice::rho - meltSnowVol * Ice::rhoSnow
            - (emp[i] - fdw[i]) * dt);

    // Effective ice salinity is always less than or equal to the SSS
    const double effectiveIceSal = std::min(Ice::s, sss[i]);
    sss[i] += +((sss[i] - effectiveIceSal) * Ice::rho * deltaIceVol // Change due to ice changes
                  + sss[i] * meltSnowVol
                  + (emp[i] - fdw[i]) * dt) // snow melt, precipitation and nudging fluxes.
        / denominator;
}

} /* namespace Nextsim */
