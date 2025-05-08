/*!
 * @file BBMDynamics.cpp
 *
 * @date 08 May 2025
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
    { BBMDynamics::CLAB_KEY, prefix + ".cohesion" },
    { BBMDynamics::SCALEC_KEY, prefix + ".scale_cohesion" },
    { BBMDynamics::REFSCALEC_KEY, prefix + ".reference_scale_cohesion" },
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
    params.compactionParam = getConfiguration(keyMap.at(C_KEY), compactionParamDefault);
    params.nu0 = getConfiguration(keyMap.at(NU_KEY), nu0Default);
    params.young = getConfiguration(keyMap.at(YOUNG_KEY), youngDefault);
    params.P0 = getConfiguration(keyMap.at(P0_KEY), P0Default);
    params.lambda0 = getConfiguration(keyMap.at(LAMBDA0_KEY), lambda0Default);
    params.alpha = getConfiguration(keyMap.at(ALPHA_KEY), alphaDefault);
    params.expPMax = getConfiguration(keyMap.at(EXPPMAX_KEY), expPMaxDefault);
    params.mu = getConfiguration(keyMap.at(MU_KEY), muDefault);
    params.comprCap = getConfiguration(keyMap.at(NMAX_KEY), comprCapDefault);
    params.cohesion = getConfiguration(keyMap.at(CLAB_KEY), cohesionDefault);
    params.scaleCohesion = getConfiguration(keyMap.at(SCALEC_KEY), scaleCohesionDefault);
    params.referenceScaleC = getConfiguration(keyMap.at(REFSCALEC_KEY), referenceScaleCDefault);
    params.nSteps = getConfiguration(keyMap.at(NSTEPS_KEY), nStepsDefault);
    params.rhoIce = getConfiguration(keyMap.at(RHOI_KEY), rhoIceDefault);
    params.rhoAtm = getConfiguration(keyMap.at(RHOA_KEY), rhoAtmDefault);
    params.rhoOcean = getConfiguration(keyMap.at(RHOO_KEY), rhoOceanDefault);
    params.CAtm = getConfiguration(keyMap.at(CATM_KEY), CAtmDefault);
    params.COcean = getConfiguration(keyMap.at(COCEAN_KEY), COceanDefault);
    params.fc = getConfiguration(keyMap.at(FC_KEY), fcDefault);
    params.oceanTurningAngle = getConfiguration(keyMap.at(ANGLE_KEY), oceanTurningAngleDefault);
}

BBMDynamics::BBMDynamics()
    : IDynamics(true)
    , kernel(params)
{
    getStore().registerArray(Protected::ICE_U, &uice, RO);
    getStore().registerArray(Protected::ICE_V, &vice, RO);

    getStore().registerArray(Protected::SHEAR, &shear, RO);
    getStore().registerArray(Protected::DIV, &divergence, RO);
    getStore().registerArray(Protected::SIGMAI, &sigmaI, RO);
    getStore().registerArray(Protected::SIGMAII, &sigmaII, RO);
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

    shear = kernel.getDG0Data(shearName);
    divergence = kernel.getDG0Data(divergenceName);
    sigmaI = kernel.getDG0Data(sigmaIName);
    sigmaII = kernel.getDG0Data(sigmaIIName);
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
            ConfigurationHelp::toString(expPMaxDefault), "",
            "Maximum compressive strength as a multiplication factor applied to the cohesion." },
        { keyMap.at(CLAB_KEY), ConfigType::NUMERIC, { "0", "∞" },
            ConfigurationHelp::toString(cohesionDefault), "Pa", "Cohesion" },
        { keyMap.at(SCALEC_KEY), ConfigType::BOOLEAN, { "true", "false" },
            ConfigurationHelp::toString(scaleCohesionDefault), "",
            "Scale the cohesion as sqrt of the grid size" },
        { keyMap.at(REFSCALEC_KEY), ConfigType::NUMERIC, { "0", "∞" },
            ConfigurationHelp::toString(referenceScaleCDefault), "m",
            "The reference scale when scaling the cohesion as sqrt of the grid size" },
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
