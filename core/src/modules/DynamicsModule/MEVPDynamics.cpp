/*!
 * @file MEVPDynamics.cpp
 *
 * @date 26 May 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Piotr Minakowski <piotr.minakowski@ovgu.de>
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#include "include/MEVPDynamics.hpp"
#include "include/constants.hpp"
#include "include/gridNames.hpp"

#include <string>
#include <vector>

namespace Nextsim {

// TODO: We should use getName() here, but it isn't static.
static const std::string prefix = "MEVPDynamics"; // MEVPDynamics::getName();
static const std::map<int, std::string> keyMap = {
    { MEVPDynamics::PSTAR_KEY, prefix + ".Pstar" },
    { MEVPDynamics::DELTA_KEY, prefix + ".DeltaMin" },
    { MEVPDynamics::C_KEY, prefix + ".C" },
    { MEVPDynamics::NSTEPS_KEY, prefix + ".nsteps" },
    { MEVPDynamics::RHOI_KEY, prefix + ".rho_ice" },
    { MEVPDynamics::RHOA_KEY, prefix + ".rho_atm" },
    { MEVPDynamics::RHOO_KEY, prefix + ".rho_ocean" },
    { MEVPDynamics::CATM_KEY, prefix + ".drag_atm" },
    { MEVPDynamics::COCEAN_KEY, prefix + ".drag_ocean" },
    { MEVPDynamics::FC_KEY, prefix + ".Coriolis_parameter" },
    { MEVPDynamics::ANGLE_KEY, prefix + ".ocean_turning_angle" },
};

void MEVPDynamics::configure()
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
}

ConfigMap MEVPDynamics::getConfiguration() const
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
    };
}

static const std::vector<std::string> namedFields = { uName, vName };
MEVPDynamics::MEVPDynamics()
    : IDynamics()
    , kernel(params)
{
}

void MEVPDynamics::setData(const ModelState::DataMap& ms)
{
    IDynamics::setData(ms);

    bool isSpherical = checkSpherical(ms);

    ModelArray coords = ms.at(coordsName);
    if (isSpherical) {
        coords *= PhysicalConstants::deg2rad;
    }
    // TODO: Some encoding of the periodic edge boundary conditions
    kernel.initialise(coords, isSpherical, ms.at(maskName));

    uice = ms.at(uName);
    vice = ms.at(vName);

    // Set the data in the kernel arrays.
    for (const auto& fieldName : namedFields) {
        kernel.setData(fieldName, ms.at(fieldName));
    }

    // Set the DG field data
    kernel.setDGArray(hiceName, hiceDG.allComponents());
    kernel.setDGArray(ciceName, ciceDG.allComponents());
}

void MEVPDynamics::update(const TimestepTime& tst)
{
    std::cout << tst.start << std::endl;

    // set the forcing velocities
    kernel.setData(uWindName, uwind);
    kernel.setData(vWindName, vwind);
    kernel.setData(uOceanName, uocean);
    kernel.setData(vOceanName, vocean);
    kernel.setData(sshName, ssh);

    kernel.update(tst);

    uice = kernel.getDG0Data(uName);
    vice = kernel.getDG0Data(vName);

    taux = kernel.getDG0Data(uIOStressName);
    tauy = kernel.getDG0Data(vIOStressName);

    shear = kernel.getDG0Data(shearName);
    divergence = kernel.getDG0Data(divergenceName);
    sigmaI = kernel.getDG0Data(sigmaIName);
    sigmaII = kernel.getDG0Data(sigmaIIName);
}

MEVPDynamics::HelpMap& MEVPDynamics::getHelpText(HelpMap& map, bool getAll)
{
    map["MEVPDynamics"] = {
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
            "Oceanic turning angle" },
    };
    return map;
}

MEVPDynamics::HelpMap& MEVPDynamics::getHelpRecursive(HelpMap& map, bool getAll)
{
    return getHelpText(map, getAll);
}

}
