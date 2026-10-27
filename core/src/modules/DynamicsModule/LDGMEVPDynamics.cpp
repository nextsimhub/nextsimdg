/*!
 * @file LDGMEVPDynamics.cpp
 *
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#include "include/LDGMEVPDynamics.hpp"
#include "include/constants.hpp"
#include "include/gridNames.hpp"

#include <stdexcept>
#include <string>
#include <vector>

namespace Nextsim {

static const std::string prefix = "LDGMEVPDynamics";
static const std::map<int, std::string> keyMap = {
    { LDGMEVPDynamics::PSTAR_KEY, prefix + ".Pstar" },
    { LDGMEVPDynamics::DELTA_KEY, prefix + ".DeltaMin" },
    { LDGMEVPDynamics::C_KEY, prefix + ".C" },
    { LDGMEVPDynamics::NSTEPS_KEY, prefix + ".nsteps" },
    { LDGMEVPDynamics::RHOI_KEY, prefix + ".rho_ice" },
    { LDGMEVPDynamics::RHOA_KEY, prefix + ".rho_atm" },
    { LDGMEVPDynamics::RHOO_KEY, prefix + ".rho_ocean" },
    { LDGMEVPDynamics::CATM_KEY, prefix + ".drag_atm" },
    { LDGMEVPDynamics::COCEAN_KEY, prefix + ".drag_ocean" },
    { LDGMEVPDynamics::FC_KEY, prefix + ".Coriolis_parameter" },
    { LDGMEVPDynamics::ANGLE_KEY, prefix + ".ocean_turning_angle" },
    { LDGMEVPDynamics::REPLACEMENT_KEY, prefix + ".replacement_pressure" },
    { LDGMEVPDynamics::ALPHA_KEY, prefix + ".alpha" },
    { LDGMEVPDynamics::BETA_KEY, prefix + ".beta" },
    { LDGMEVPDynamics::LDGA_KEY, prefix + ".ldg_a" },
    { LDGMEVPDynamics::LDGB_KEY, prefix + ".ldg_b" },
    { LDGMEVPDynamics::PCGTOL_KEY, prefix + ".pcg_tolerance" },
    { LDGMEVPDynamics::PCGMAXIT_KEY, prefix + ".pcg_max_iterations" },
};

static const LDGParameters ldgDefaults;

void LDGMEVPDynamics::configure()
{
    Module::Module<Nextsim::IDamageHealing>::setImplementation("Nextsim::NoHealing");

    params.pStar = Configured::getConfiguration(keyMap.at(PSTAR_KEY), pStarDefault);
    params.deltaMin = Configured::getConfiguration(keyMap.at(DELTA_KEY), deltaMinDefault);
    params.compactionParam = Configured::getConfiguration(keyMap.at(C_KEY), compactionParamDefault);
    params.nSteps = Configured::getConfiguration(keyMap.at(NSTEPS_KEY), nStepsDefault);
    params.rhoIce = Configured::getConfiguration(keyMap.at(RHOI_KEY), rhoIceDefault);
    params.rhoAtm = Configured::getConfiguration(keyMap.at(RHOA_KEY), rhoAtmDefault);
    params.rhoOcean = Configured::getConfiguration(keyMap.at(RHOO_KEY), rhoOceanDefault);
    params.CAtm = Configured::getConfiguration(keyMap.at(CATM_KEY), CAtmDefault);
    params.COcean = Configured::getConfiguration(keyMap.at(COCEAN_KEY), COceanDefault);
    params.fc = Configured::getConfiguration(keyMap.at(FC_KEY), fcDefault);
    params.oceanTurningAngle
        = Configured::getConfiguration(keyMap.at(ANGLE_KEY), oceanTurningAngleDefault);
    params.replacementPressure = Configured::getConfiguration(keyMap.at(REPLACEMENT_KEY), false);

    ldgParams.alpha = Configured::getConfiguration(keyMap.at(ALPHA_KEY), ldgDefaults.alpha);
    ldgParams.beta = Configured::getConfiguration(keyMap.at(BETA_KEY), ldgDefaults.beta);
    ldgParams.a = Configured::getConfiguration(keyMap.at(LDGA_KEY), ldgDefaults.a);
    ldgParams.b = Configured::getConfiguration(keyMap.at(LDGB_KEY), ldgDefaults.b);
    ldgParams.pcgTolerance
        = Configured::getConfiguration(keyMap.at(PCGTOL_KEY), ldgDefaults.pcgTolerance);
    ldgParams.pcgMaxIterations
        = Configured::getConfiguration(keyMap.at(PCGMAXIT_KEY), ldgDefaults.pcgMaxIterations);
}

ConfigMap LDGMEVPDynamics::getConfiguration() const
{
    return {
        { keyMap.at(PSTAR_KEY), params.pStar },
        { keyMap.at(DELTA_KEY), params.deltaMin },
        { keyMap.at(C_KEY), params.compactionParam },
        { keyMap.at(NSTEPS_KEY), params.nSteps },
        { keyMap.at(RHOI_KEY), params.rhoIce },
        { keyMap.at(RHOA_KEY), params.rhoAtm },
        { keyMap.at(RHOO_KEY), params.rhoOcean },
        { keyMap.at(CATM_KEY), params.CAtm },
        { keyMap.at(COCEAN_KEY), params.COcean },
        { keyMap.at(FC_KEY), params.fc },
        { keyMap.at(ANGLE_KEY), params.oceanTurningAngle },
        { keyMap.at(REPLACEMENT_KEY), params.replacementPressure },
        { keyMap.at(ALPHA_KEY), ldgParams.alpha },
        { keyMap.at(BETA_KEY), ldgParams.beta },
        { keyMap.at(LDGA_KEY), ldgParams.a },
        { keyMap.at(LDGB_KEY), ldgParams.b },
        { keyMap.at(PCGTOL_KEY), ldgParams.pcgTolerance },
        { keyMap.at(PCGMAXIT_KEY), ldgParams.pcgMaxIterations },
    };
}

static const std::vector<std::string> namedFields = { uName, vName };
LDGMEVPDynamics::LDGMEVPDynamics()
    : IDynamics()
    , kernel(params, ldgParams)
{
}

void LDGMEVPDynamics::prepareAdvection() { kernel.prepareAdvection(); }

void LDGMEVPDynamics::setData(const ModelState::DataMap& ms)
{
    IDynamics::setData(ms);

    const bool isSpherical = checkSpherical(ms);

    ModelArray coords = ms.at(coordsName);
    if (isSpherical) {
        coords *= PhysicalConstants::deg2rad;
    }

    kernel.initialise(coords, isSpherical, ms.at(maskName));

    uiceAccessor.getHostRW() = ms.at(uName);
    viceAccessor.getHostRW() = ms.at(vName);

    // Set the data in the kernel arrays.
    for (const auto& fieldName : namedFields) {
        kernel.setData(fieldName, ms.at(fieldName));
    }
}

void LDGMEVPDynamics::update(const TimestepTime& tst)
{
    std::cout << tst.start << std::endl;

    // set dg fields
    kernel.setDGArray(hiceName, hiceDGAccessor.getAutoRW());
    kernel.setDGArray(ciceName, ciceDGAccessor.getAutoRW());
    kernel.setDGArray(hsnowName, hsnowDGAccessor.getAutoRW());

    // set the forcing velocities
    kernel.setData(uWindName, uwindAccessor.getAutoRO());
    kernel.setData(vWindName, vwindAccessor.getAutoRO());
    kernel.setData(uOceanName, uoceanAccessor.getAutoRO());
    kernel.setData(vOceanName, voceanAccessor.getAutoRO());
    kernel.setData(sshName, sshAccessor.getHostRO());

    kernel.update(tst);

    uiceAccessor.getHostRW() = kernel.getDG0Data(uName);
    viceAccessor.getHostRW() = kernel.getDG0Data(vName);

    tauxAccessor.getHostRW() = kernel.getDG0Data(uIOStressName);
    tauyAccessor.getHostRW() = kernel.getDG0Data(vIOStressName);

    shearAccessor.getHostRW() = kernel.getDG0Data(shearName);
    divergenceAccessor.getHostRW() = kernel.getDG0Data(divergenceName);
    sigmaIAccessor.getHostRW() = kernel.getDG0Data(sigmaIName);
    sigmaIIAccessor.getHostRW() = kernel.getDG0Data(sigmaIIName);
}

void LDGMEVPDynamics::advectField(
    FloatType timestep, ModelArray& field, FloatType lowerLimit, FloatType upperLimit)
{
    kernel.advectField(timestep, field, lowerLimit, upperLimit);
}

#ifdef USE_KOKKOS
void LDGMEVPDynamics::advectField(
    FloatType timestep, const DeviceViewMA& field, FloatType lowerLimit, FloatType upperLimit)
{
    throw std::runtime_error("LDGMEVPDynamics: Kokkos is not supported");
}
#endif

LDGMEVPDynamics::HelpMap& LDGMEVPDynamics::getHelpText(HelpMap& map, bool getAll)
{
    map["LDGMEVPDynamics"] = {
        { keyMap.at(PSTAR_KEY), ConfigType::NUMERIC, { "0", "∞" },
            ConfigurationHelp::toString(pStarDefault), "Pa", "The sea-ice strength parameter P*" },
        { keyMap.at(DELTA_KEY), ConfigType::NUMERIC, { "0", "∞" },
            ConfigurationHelp::toString(deltaMinDefault), "1/s",
            "Capping value of the 𝛥 variable" },
        { keyMap.at(C_KEY), ConfigType::NUMERIC, { "-∞", "0" },
            ConfigurationHelp::toString(compactionParamDefault), "[None]",
            "The compaction parameter C" },
        { keyMap.at(NSTEPS_KEY), ConfigType::NUMERIC, { "1", "∞" },
            ConfigurationHelp::toString(nStepsDefault), "[No unit]",
            "The number of sub-cycling steps" },
        { keyMap.at(RHOI_KEY), ConfigType::NUMERIC, { "0", "∞" },
            ConfigurationHelp::toString(rhoIceDefault), "kg/m^3", "Density of sea ice" },
        { keyMap.at(RHOA_KEY), ConfigType::NUMERIC, { "0", "∞" },
            ConfigurationHelp::toString(rhoAtmDefault), "kg/m^3", "Density of air" },
        { keyMap.at(RHOO_KEY), ConfigType::NUMERIC, { "0", "∞" },
            ConfigurationHelp::toString(rhoOceanDefault), "kg/m^3", "Density of ocean" },
        { keyMap.at(CATM_KEY), ConfigType::NUMERIC, { "0", "∞" },
            ConfigurationHelp::toString(CAtmDefault), "[No unit]",
            "Ice-atmosphere drag coefficient" },
        { keyMap.at(COCEAN_KEY), ConfigType::NUMERIC, { "0", "∞" },
            ConfigurationHelp::toString(COceanDefault), "[No unit]", "Ice-ocean drag coefficient" },
        { keyMap.at(FC_KEY), ConfigType::NUMERIC, { "0", "∞" },
            ConfigurationHelp::toString(fcDefault), "[No unit]",
            "Coriolis parameter (constant across the domain)" },
        { keyMap.at(ANGLE_KEY), ConfigType::NUMERIC, { "0", "90" },
            ConfigurationHelp::toString(oceanTurningAngleDefault), "degrees",
            "Oceanic turning angle (ice-ocean stress output only)" },
        { keyMap.at(REPLACEMENT_KEY), ConfigType::BOOLEAN, { "true", "false" }, "false", "",
            "Use the replacement pressure P*Delta/(Delta+DeltaMin) in the stress" },
        { keyMap.at(ALPHA_KEY), ConfigType::NUMERIC, { "1", "∞" },
            ConfigurationHelp::toString(ldgDefaults.alpha), "[No unit]",
            "mEVP relaxation parameter of the stress" },
        { keyMap.at(BETA_KEY), ConfigType::NUMERIC, { "0", "∞" },
            ConfigurationHelp::toString(ldgDefaults.beta), "[No unit]",
            "mEVP relaxation parameter of the velocity" },
        { keyMap.at(LDGA_KEY), ConfigType::NUMERIC, { "-0.5", "0.5" },
            ConfigurationHelp::toString(ldgDefaults.a), "[No unit]",
            "LDG alternating flux: u_hat = (0.5+a) u1 + (0.5-a) u2" },
        { keyMap.at(LDGB_KEY), ConfigType::NUMERIC, { "0", "∞" },
            ConfigurationHelp::toString(ldgDefaults.b), "Pa s",
            "LDG penalty of velocity jumps, scaled with 1/h" },
        { keyMap.at(PCGTOL_KEY), ConfigType::NUMERIC, { "0", "1" },
            ConfigurationHelp::toString(ldgDefaults.pcgTolerance), "[No unit]",
            "Relative residual of the PCG solver for the velocity" },
        { keyMap.at(PCGMAXIT_KEY), ConfigType::INTEGER, { "1", "∞" },
            ConfigurationHelp::toString(ldgDefaults.pcgMaxIterations), "[No unit]",
            "Maximum number of PCG iterations for the velocity" },
    };
    return map;
}

LDGMEVPDynamics::HelpMap& LDGMEVPDynamics::getHelpRecursive(HelpMap& map, bool getAll)
{
    return getHelpText(map, getAll);
}

}
