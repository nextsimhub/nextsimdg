/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#include "include/SlabOcean.hpp"

#include "include/constants.hpp"
#include "include/gridNames.hpp"

#include "include/KernelAlternatives.hpp"
#include "kokkos/include/KokkosTimer.hpp"

#include <map>
#include <string>

namespace Nextsim {

const FloatType SlabOcean::defaultRelaxationTime = 30; // 30 days in seconds

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
            ConfigurationHelp::toString(defaultRelaxationTime), "days",
            "Relaxation time of the slab ocean to external temperature forcing." },
        { keyMap.at(TIMES_KEY), ConfigType::NUMERIC, { "0", "∞" },
            ConfigurationHelp::toString(defaultRelaxationTime), "days",
            "Relaxation time of the slab ocean to external salinity forcing." },
    };
    return map;
};

void SlabOcean::setData(const ModelState::DataMap& ms)
{
    HField& qdw = qdwAccessor.getHostRW();
    qdw.reinitialize();
    HField& fdw = fdwAccessor.getHostRW();
    fdw.reinitialize();
    HField& sstSlab = sstSlabAccessor.getHostRW();
    sstSlab.reinitialize();
    HField& sssSlab = sssSlabAccessor.getHostRW();
    sssSlab.reinitialize();
}

void SlabOcean::update(const TimestepTime& tst)
{
    static KokkosTimer<true> timer("SlabOcean");

    timer.start();

    auto execSpace = DefaultExecutionSpace();
    auto& sstSlab = sstSlabAccessor.getAutoRW(execSpace);
    auto& fdw = fdwAccessor.getAutoRW(execSpace);
    auto& qdw = qdwAccessor.getAutoRW(execSpace);
    auto& sssSlab = sssSlabAccessor.getAutoRW(execSpace);
    const auto& fwFlux = fwFluxAccessor.getAutoRO(execSpace);
    const auto& sFlux = sFluxAccessor.getAutoRO(execSpace);
    const auto& sssExt = sssExtAccessor.getAutoRO(execSpace);
    const auto& sst = sstAccessor.getAutoRO(execSpace);
    const auto& sstExt = sstExtAccessor.getAutoRO(execSpace);
    const auto& qswNet = qswNetAccessor.getAutoRO(execSpace);
    const auto& sss = sssAccessor.getAutoRO(execSpace);
    const auto& qNoSun = qNoSunAccessor.getAutoRO(execSpace);
    const auto& cpml = cpmlAccessor.getAutoRO(execSpace);

    const FloatType dt = tst.step.seconds();
    const FloatType rRelaxationTimeT = 1 / (relaxationTimeT * 86400);
    const FloatType rRelaxationTimeS = 1 / (relaxationTimeS * 86400);

    overElementsAuto(OVER_ELEMENTS_LAMBDA(const ElementIndex i) {
        // Slab SST update
        qdw[i] = (sstExt[i] - sst[i]) * cpml[i] * rRelaxationTimeT;
        sstSlab[i] = sst[i] - dt * (qswNet[i] + qNoSun[i] - qdw[i]) / cpml[i];

        // Slab SSS update
        const FloatType arealDensity
            = cpml[i] / Water::cp; // density times depth, or cpml divided by cp
        /* Just use a salt flux as the nudging flux. This is simplified compared to the
         * finiteelement.cpp calculation
         * Fdw = delS * mld * physical::rhow /(timeS*M_sss[i] - ddt*delS)
         * where delS = sssSlab - sssExt
         */
        fdw[i] = (sssExt[i] - sss[i]) * arealDensity * rRelaxationTimeS;

        // Mass per unit area after all the changes in water volume
        // sFlux is in kg/m^2/s, but we need PSU/m^2/s
        sssSlab[i] = (sss[i] * arealDensity * (fdw[i] - 1e3_ft * sFlux[i]) * dt)
            / (arealDensity - fwFlux[i] * dt);
    });
    timer.stop();
}

} /* namespace Nextsim */
