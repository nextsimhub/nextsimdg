/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#include "include/SlabOcean.hpp"

#include "include/constants.hpp"
#include "include/gridNames.hpp"

#include <map>
#include <string>

namespace Nextsim {

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
}

ConfigMap SlabOcean::getConfiguration() const
{
    return {
        { keyMap.at(TIMET_KEY), relaxationTimeT },
        { keyMap.at(TIMES_KEY), relaxationTimeS },
    };
}

ModelState SlabOcean::getStatePrognostic() const
{
    return { {
                 { sstName, sstSlabAccessor.getHostRO() },
                 { sssName, sssSlabAccessor.getHostRO() },
             },
        getConfiguration() };
}

ModelState SlabOcean::getStateDiagnostic() const
{
    ModelState state = { {
                             { "Q_slab", qdwAccessor.getHostRO() },
                             { "F_slab", fdwAccessor.getHostRO() },
                         },
        {} };

    return state.merge(getStatePrognostic());
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
    HField& qdw = qdwAccessor.getHostRW();
    qdw.resize();
    HField& fdw = fdwAccessor.getHostRW();
    fdw.resize();
    HField& sstSlab = sstSlabAccessor.getHostRW();
    sstSlab.resize();
    HField& sssSlab = sssSlabAccessor.getHostRW();
    sssSlab.resize();
}

void SlabOcean::update(const TimestepTime& tst)
{
    HField& sstSlab = sstSlabAccessor.getHostRW();
    HField& fdw = fdwAccessor.getHostRW();
    HField& qdw = qdwAccessor.getHostRW();
    HField& sssSlab = sssSlabAccessor.getHostRW();
    const HField& fwFlux = fwFluxAccessor.getHostRO();
    const HField& sssExt = sssExtAccessor.getHostRO();
    const HField& sst = sstAccessor.getHostRO();
    const HField& sstExt = sstExtAccessor.getHostRO();
    const HField& qswNet = qswNetAccessor.getHostRO();
    const HField& sss = sssAccessor.getHostRO();
    const HField& qNoSun = qNoSunAccessor.getHostRO();
    const HField& cpml = cpmlAccessor.getHostRO();

    dt = tst.step.seconds();
    overElements(
        [&](const size_t i, const TimestepTime& tsTime) {
            // Slab SST update
            qdw[i] = (sstExt[i] - sst[i]) * cpml[i] / relaxationTimeT;
            sstSlab[i] = sst[i] - dt * (qswNet[i] + qNoSun[i] - qdw[i]) / cpml[i];

            // Slab SSS update
            const double arealDensity
                = cpml[i] / Water::cp; // density times depth, or cpml divided by cp
            // This is simplified compared to the finiteelement.cpp calculation
            // Fdw = delS * mld * physical::rhow /(timeS*M_sss[i] - ddt*delS) where delS = sssSlab -
            // sssExt
            fdw[i] = (1 - sssExt[i] / sss[i]) * arealDensity / relaxationTimeS;

            // Mass per unit area after all the changes in water volume
            // Clamp the denominator to be at least 1 m deep, i.e. at least Water::rho kg m⁻²
            const double denominator
                = std::max(arealDensity - (fwFlux[i] - fdw[i]) * dt, Water::rhoOcean);
            sssSlab[i] = sss[i] + (sss[i] * fwFlux[i] - fdw[i] * dt) / denominator;
        },
        tst);
}
/*
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
}*/

} /* namespace Nextsim */
