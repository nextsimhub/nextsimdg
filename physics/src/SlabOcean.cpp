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

const FloatType SlabOcean::defaultRelaxationTime = 7; // Unit is days

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

    // Compute in double because the flux computations are sensitive to catastrophic cancellation.
    const double dt = tst.step.seconds();
    const double rRelaxationTimeT = 1.0 / (relaxationTimeT * 86400.0);
    const double rRelaxationTimeS = 1.0 / (relaxationTimeS * 86400.0);

    overElementsAuto(OVER_ELEMENTS_LAMBDA(const ElementIndex i) {
        // Slab SST update
        const double sstExtD = static_cast<double>(sstExt[i]);
        const double sstD = static_cast<double>(sst[i]);
        const double cpmlD = static_cast<double>(cpml[i]);
        const double qswNetD = static_cast<double>(qswNet[i]);
        const double qNoSunD = static_cast<double>(qNoSun[i]);

        const double qdwD = (sstExtD - sstD) * cpmlD * rRelaxationTimeT;
        qdw[i] = qdwD;
        sstSlab[i] = sstD - dt * (qswNetD + qNoSunD - qdwD) / cpmlD;

        // Slab SSS update
        const double sssD = static_cast<double>(sss[i]);
        const double sssExtD = static_cast<double>(sssExt[i]);
        const double sFluxD = static_cast<double>(sFlux[i]);
        const double fwFluxD = static_cast<double>(fwFlux[i]);

        const double arealDensity = cpmlD / Water::cp; // density times depth, or cpml divided by cp
        /* Just use a salt flux as the nudging flux. This is simplified compared to the
         * finiteelement.cpp calculation
         * Fdw = delS * mld * physical::rhow /(timeS*M_sss[i] - ddt*delS)
         * where delS = sssSlab - sssExt
         */
        const double fdwD = (sssExtD - sssD) * arealDensity * rRelaxationTimeS;
        fdw[i] = fdwD;

        // Mass per unit area after all the changes in water volume
        // sFlux is in kg/m^2/s, but we need PSU/m^2/s
        sssSlab[i]
            = (sssD * arealDensity + (fdwD - 1e3 * sFluxD) * dt) / (arealDensity - fwFluxD * dt);
    });
    timer.stop();
}

} /* namespace Nextsim */
