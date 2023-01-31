/*!
 * @file SlabOcean.cpp
 *
 * @date 27 Jan 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/SlabOcean.hpp"

#include "include/constants.hpp"

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

    // Mirror the model SS* data into the slab SS* data
    registerProtectedArray(ProtectedArray::SLAB_SST, &sstSlab.data());
    registerProtectedArray(ProtectedArray::SLAB_SST, &sstSlab.data());

    registerProtectedArray(ProtectedArray::SLAB_QDW, &qdw);
    registerProtectedArray(ProtectedArray::SLAB_FDW, &fdw);
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
    qdw = (sstExt - sstSlab) * cpml / timeT;
    sstSlab -= dt * (qio + qow -qdw) / cpml;
    // Slab SSS update
    HField arealDensity = cpml / Water::cp; // density times depth, or cpml divided by cp
    // This is simplified compared to the finiteelement.cpp calculation
    // Fdw = delS * mld * physical::rhow /(timeS*M_sss[i] - ddt*delS) where delS = sssSlab - sssExt
    fdw = ( 1- sssExt / sssSlab) * arealDensity / timeS;
    // ice volume change, both laterally and vertically
    HField deltaIceVol = newIce + deltaHice * cice;
    // change in snow volume due to melting (should be < 0)
    HField meltSnowVol = deltaSmelt * cice;
    // Mass per unit area after all the changes in water volume
    HField denominator = arealDensity - deltaIceVol * Ice::rho - meltSnowVol * Ice::rhoSnow - (emp - fdw) * dt;
    // Clamp the denominator to be at least 1 m deep, i.e. at least Water::rho kg m⁻²
    denominator.clampAbove(Water::rho);
    // Effective ice salinity is always less than or equal to the SSS
    HField effectiveIceSal = sssSlab;
    effectiveIceSal.clampBelow(Ice::s);
    sssSlab += ( (sssSlab - effectiveIceSal) * Ice::rho * deltaIceVol // Change due to ice changes
            + sssSlab * meltSnowVol + (emp - fdw) * dt) // snow melt, precipitation and nudging fluxes.
                    / denominator;
}

} /* namespace Nextsim */
