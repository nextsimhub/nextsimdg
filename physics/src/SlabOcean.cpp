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
    static KokkosTimer<true> timer("SlabOcean");

    timer.start();

    auto execSpace = DefaultExecutionSpace();
    auto& sstSlab = sstSlabAccessor.getAutoRW(execSpace);
    auto& fdw = fdwAccessor.getAutoRW(execSpace);
    auto& qdw = qdwAccessor.getAutoRW(execSpace);
    auto& sssSlab = sssSlabAccessor.getAutoRW(execSpace);
    const auto& fwFlux = fwFluxAccessor.getAutoRO(execSpace);
    const auto& sssExt = sssExtAccessor.getAutoRO(execSpace);
    const auto& sst = sstAccessor.getAutoRO(execSpace);
    const auto& sstExt = sstExtAccessor.getAutoRO(execSpace);
    const auto& qswNet = qswNetAccessor.getAutoRO(execSpace);
    const auto& sss = sssAccessor.getAutoRO(execSpace);
    const auto& qNoSun = qNoSunAccessor.getAutoRO(execSpace);
    const auto& cpml = cpmlAccessor.getAutoRO(execSpace);

    const double dt = tst.step.seconds();
    const double relaxationTimeT = SlabOcean::relaxationTimeT;
    const double relaxationTimeS = SlabOcean::relaxationTimeS;
    
    overElementsAuto(OVER_ELEMENTS_LAMBDA(const ElementIndex i) {
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

        // the device compiler does not like a global constant appearing in the argument list of
        // a template function: "Water::rhoOcean" is undefined in device code"
        const double rhoOcean = Water::rhoOcean;
        // Mass per unit area after all the changes in water volume
        // Clamp the denominator to be at least 1 m deep, i.e. at least Water::rho kg m⁻²
        const double denominator = Utils::max(arealDensity - (fwFlux[i] - fdw[i]) * dt, rhoOcean);
        sssSlab[i] = sss[i] + (sss[i] * fwFlux[i] - fdw[i] * dt) / denominator;
    });
    timer.stop();
}

} /* namespace Nextsim */
