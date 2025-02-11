/*!
 * @file BBMDynamics.cpp
 *
 * @date 05 Dec 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#include "include/BBMDynamics.hpp"
#include "include/constants.hpp"
#include "include/gridNames.hpp"

namespace Nextsim {

static const std::vector<std::string> namedFields = { hiceName, ciceName, uName, vName };
static const std::map<std::string, std::pair<ModelArray::Type, double>> defaultFields = {
    { damageName, { ModelArray::Type::H, 1.0 } },
};

// TODO: We should use getName() here, but it isn't static.
static const std::string prefix = "BBMDynamics"; // MEVPDynamics::getName();
static const std::map<int, std::string> keyMap = {
    { BBMDynamics::C_KEY, prefix + ".C" },
    { BBMDynamics::NU_KEY, prefix + ".nu" },
    { BBMDynamics::YOUNG_KEY, prefix + ".young" },
    { BBMDynamics::P0_KEY, prefix + ".P0" },
    { BBMDynamics::LAMBDA0_KEY, prefix + ".lambda0" },
    { BBMDynamics::ALPHA_KEY, prefix + ".alpha" },
    { BBMDynamics::EXPPMAX_KEY, prefix + ".exppmax" },
    { BBMDynamics::MU_KEY, prefix + ".mu" },
    { BBMDynamics::NMAX_KEY, prefix + ".namx" },
    { BBMDynamics::CLAB_KEY, prefix + ".clab" },
    { BBMDynamics::NSTEPS_KEY, prefix + ".nsteps" },
    { BBMDynamics::RHOI_KEY, prefix + ".rho_ice" },
    { BBMDynamics::RHOA_KEY, prefix + ".rho_atm" },
    { BBMDynamics::RHOO_KEY, prefix + ".rho_ocean" },
    { BBMDynamics::CATM_KEY, prefix + ".drag_atm" },
    { BBMDynamics::COCEAN_KEY, prefix + ".drag_ocean" },
    { BBMDynamics::FC_KEY, prefix + ".Coriolis_parameter" },
    { BBMDynamics::ANGLE_KEY, prefix + ".ocean_turning_angle" },
};

void BBMDynamics::configure()
{
    params.compactionParam = Configured::getConfiguration(keyMap.at(C_KEY), compactionParamDefault);
    params.nu0 = Configured::getConfiguration(keyMap.at(NU_KEY), nu0Default);
    params.young = Configured::getConfiguration(keyMap.at(YOUNG_KEY), youngDefault);
    params.P0 = Configured::getConfiguration(keyMap.at(P0_KEY), P0Default);
    params.lambda0 = Configured::getConfiguration(keyMap.at(LAMBDA0_KEY), lambda0Default);
    params.alpha = Configured::getConfiguration(keyMap.at(ALPHA_KEY), alphaDefault);
    params.expPMax = Configured::getConfiguration(keyMap.at(EXPPMAX_KEY), expPMaxDefault);
    params.mu = Configured::getConfiguration(keyMap.at(MU_KEY), muDefault);
    params.comprCap = Configured::getConfiguration(keyMap.at(NMAX_KEY), comprCapDefault);
    params.cLab = Configured::getConfiguration(keyMap.at(CLAB_KEY), cLabDefault);
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

BBMDynamics::BBMDynamics()
    : IDynamics(true)
    , kernel(params)
{
    getStore().registerArray(Protected::ICE_U, &uice, RO);
    getStore().registerArray(Protected::ICE_V, &vice, RO);
}

void BBMDynamics::setData(const ModelState::DataMap& ms)
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
    // Required data
    for (const auto& fieldName : namedFields) {
        kernel.setData(fieldName, ms.at(fieldName));
    }
    // Data that can have a default value
    for (const auto entry : defaultFields) {
        // Directly add data that is supplied
        const std::string& fieldName = entry.first;
        if (ms.count(fieldName) > 0) {
            kernel.setData(fieldName, ms.at(fieldName));
        } else {
            // Fill data that is not supplied, masking if the mask is available
            ModelArray data(entry.second.first);
            data.resize();
            // Fill the default value
            data = entry.second.second;
            // Mask the default data
            kernel.setData(fieldName, mask(data));
        }
    }
}

void BBMDynamics::update(const TimestepTime& tst)
{
    std::cout << tst.start << std::endl;

    // Fill the updated damage array with the initial value
    damage = damage0.data();

    // set the updated ice thickness, concentration and damage
    kernel.setData(hiceName, hice.data());
    kernel.setData(ciceName, cice.data());
    kernel.setData(damageName, damage);

    // set the forcing velocities
    kernel.setData(uWindName, uwind.data());
    kernel.setData(vWindName, vwind.data());
    kernel.setData(uOceanName, uocean.data());
    kernel.setData(vOceanName, vocean.data());
    kernel.setData(sshName, ssh.data());

    /*
     * Ice velocity components are stored in the dynamics, and not changed by the model outside the
     * dynamics kernel. Hence they are not set at this point.
     */

    kernel.update(tst);

    hice.data() = kernel.getDG0Data(hiceName);
    cice.data() = kernel.getDG0Data(ciceName);
    damage = kernel.getDG0Data(damageName);

    uice = kernel.getDG0Data(uName);
    vice = kernel.getDG0Data(vName);

    taux = kernel.getDG0Data(uIOStressName);
    tauy = kernel.getDG0Data(vIOStressName);
}

// All data for prognostic output
ModelState BBMDynamics::getState() const
{
    // Get the velocities from IDynamics
    ModelState state(IDynamics::getState());

    // Kernel prognostic fields
    state.merge({
        { hiceName, kernel.getDGData(hiceName) },
        { ciceName, kernel.getDGData(ciceName) },
        { damageName, kernel.getDGData(damageName) },
    });

    return state;
}

ModelState BBMDynamics::getStateRecursive(const OutputSpec& os) const
{
    // Base class state
    ModelState state(IDynamics::getStateRecursive(os));

    if (os.allComponents()) {
        state.merge({
            { hiceName, kernel.getDGData(hiceName) },
            { ciceName, kernel.getDGData(ciceName) },
            { damageName, kernel.getDGData(damageName) },
        });
    }
    return state;
}

BBMDynamics::HelpMap& BBMDynamics::getHelpText(HelpMap& map, bool getAll)
{
    map["BBMDynamics"] = {
        { keyMap.at(C_KEY), ConfigType::NUMERIC, { "-∞", "0" },
            ConfigurationHelp::toString(compactionParamDefault), "[None]",
            "The compaction parameter C" },
        { keyMap.at(NU_KEY), ConfigType::NUMERIC, { "-∞", "0" },
            ConfigurationHelp::toString(nu0Default), "[None]", "Poisson's ratio, 𝜈" },
        { keyMap.at(YOUNG_KEY), ConfigType::NUMERIC, { "0", "∞" },
            ConfigurationHelp::toString(youngDefault), "Pa", "Young's modulus, Y" },
        { keyMap.at(P0_KEY), ConfigType::NUMERIC, { "0", "∞" },
            ConfigurationHelp::toString(P0Default), "Pa", "Ice strength scaling parameter" },
        { keyMap.at(LAMBDA0_KEY), ConfigType::NUMERIC, { "0", "∞" },
            ConfigurationHelp::toString(lambda0Default), "s", "Undamaged relaxation time scale" },
        { keyMap.at(ALPHA_KEY), ConfigType::NUMERIC, { "2", "∞" },
            ConfigurationHelp::toString(alphaDefault), "[None]", "Damage parameter" },
        { keyMap.at(EXPPMAX_KEY), ConfigType::NUMERIC, { "0", "2" },
            ConfigurationHelp::toString(expPMaxDefault), "[None]",
            "Exponent for thickness scaling of P_{max}" },
        { keyMap.at(MU_KEY), ConfigType::NUMERIC, { "0", "1" },
            ConfigurationHelp::toString(expPMaxDefault), "[None]",
            "Internal friction coefficient, 𝜇" },
        { keyMap.at(NMAX_KEY), ConfigType::NUMERIC, { "0", "∞" },
            ConfigurationHelp::toString(expPMaxDefault), "Pa",
            "Maximum compressive strength (at the lab scale)" },
        { keyMap.at(CLAB_KEY), ConfigType::NUMERIC, { "0", "∞" },
            ConfigurationHelp::toString(cLabDefault), "Pa", "Cohesion (at the lab scale)" },
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

BBMDynamics::HelpMap& BBMDynamics::getHelpRecursive(HelpMap& map, bool getAll)
{
    return getHelpText(map, getAll);
}

} /* namespace Nextsim */
