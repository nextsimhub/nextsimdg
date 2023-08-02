/*!
 * @file IceGrowth.cpp
 *
 * @date Mar 15, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/IceGrowth.hpp"

#include "include/Module.hpp"
#include "include/constants.hpp"

namespace Nextsim {

template <>
const std::map<int, std::string> Configured<IceGrowth>::keyMap = {
    { IceGrowth::ICE_THERMODYNAMICS_KEY, "IceThermodynamicsModel" },
    { IceGrowth::LATERAL_GROWTH_KEY, "LateralIceModel" },
    { IceGrowth::MINC_KEY, "nextsim_thermo.min_conc" },
    { IceGrowth::MINH_KEY, "nextsim_thermo.min_thick" },
};

IceGrowth::IceGrowth()
    : hice(ModelArray::Type::H)
    , cice(ModelArray::Type::H)
    , hsnow(ModelArray::Type::H)
    , hice0(ModelArray::Type::H)
    , hsnow0(ModelArray::Type::H)
    , newice(ModelArray::Type::H)
    , deltaCIce(ModelArray::Type::H)
    , deltaCFreeze(ModelArray::Type::H)
    , deltaCMelt(ModelArray::Type::H)
    , hIceCell(getProtectedArray())
    , hSnowCell(getProtectedArray())
    , cice0(getProtectedArray())
    , qow(getSharedArray())
    , mixedLayerBulkHeatCapacity(getProtectedArray())
    , sst(getProtectedArray())
    , tf(getProtectedArray())
    , deltaHi(getSharedArray())
{
    registerModule();
    registerSharedArray(SharedArray::H_ICE, &hice);
    registerSharedArray(SharedArray::C_ICE, &cice);
    registerSharedArray(SharedArray::H_SNOW, &hsnow);
    registerSharedArray(SharedArray::NEW_ICE, &newice);
    registerSharedArray(SharedArray::HSNOW_MELT, &snowMelt);
    registerSharedArray(SharedArray::DELTA_CICE, &deltaCIce);

    registerProtectedArray(ProtectedArray::HTRUE_ICE, &hice0);
    registerProtectedArray(ProtectedArray::HTRUE_SNOW, &hsnow0);
}

void IceGrowth::setData(const ModelState::DataMap& ms)
{
    iVertical->setData(ms);
    iLateral->setData(ms);

    hice.resize();
    cice.resize();
    hsnow.resize();
    hice0.resize();
    hsnow0.resize();
    newice.resize();
    snowMelt.resize();
    deltaCFreeze.resize();
    deltaCMelt.resize();
    deltaCIce.resize();
}

ModelState IceGrowth::getState() const
{
    return { {
                 { "hice_updated", hice },
                 { "cice_updated", cice },
                 { "hsnow_updated", hsnow },
                 { "hice_initial", hice0 },
                 { "cice_initial", cice0 },
                 { "hsnow_initial", hsnow0 },
             },
        getConfiguration() };
}

ModelState IceGrowth::getStateRecursive(const OutputSpec& os) const
{
    ModelState state(getState());
    // Merge in other states here
    state.merge(iLateral->getStateRecursive(os));
    state.merge(iVertical->getStateRecursive(os));

    return os ? state : ModelState();
}

IceGrowth::HelpMap& IceGrowth::getHelpText(HelpMap& map, bool getAll)
{
    map["IceGrowth"] = {
        { keyMap.at(MINC_KEY), ConfigType::NUMERIC, { "0", "1" },
            std::to_string(IceMinima::cMinDefault), "", "Minimum allowed ice concentration." },
        { keyMap.at(MINH_KEY), ConfigType::NUMERIC, { "0", "∞" },
            std::to_string(IceMinima::hMinDefault), "m", "Minimum allowed ice thickness." },
    };
    return map;
}
IceGrowth::HelpMap& IceGrowth::getHelpRecursive(HelpMap& map, bool getAll)
{
    getHelpText(map, getAll);
    Module::getHelpRecursive<IIceThermodynamics>(map, getAll);
    Module::getHelpRecursive<ILateralIceSpread>(map, getAll);
    return map;
}

void IceGrowth::configure()
{
    // Configure constants
    IceMinima::cMin = Configured::getConfiguration(keyMap.at(MINC_KEY), IceMinima::cMinDefault);
    IceMinima::hMin = Configured::getConfiguration(keyMap.at(MINH_KEY), IceMinima::hMinDefault);

    // Configure the vertical and lateral growth modules
    iVertical = std::move(Module::getInstance<IIceThermodynamics>());
    iLateral = std::move(Module::getInstance<ILateralIceSpread>());
    tryConfigure(*iVertical);
    tryConfigure(*iLateral);
}

ConfigMap IceGrowth::getConfiguration() const
{
    return {
        { keyMap.at(MINC_KEY), IceMinima::cMin },
        { keyMap.at(MINH_KEY), IceMinima::hMin },
    };
}

void IceGrowth::update(const TimestepTime& tsTime)
{

    // Copy the ice data from the prognostic fields to the modifiable fields.
    initializeThicknesses();
    overElements(
            std::bind(&IceGrowth::applyLimits, this, std::placeholders::_1, std::placeholders::_2),
            tsTime);
    std::cerr << "After thick init: hice = " << hice(79,67) << " cice = " << cice(79,67) << std::endl;

    iVertical->update(tsTime);
    std::cerr << "After thermo: hice = " << hice(79,67) << " cice = " << cice(79,67) << std::endl;
    // new ice formation
    overElements(
        std::bind(&IceGrowth::updateWrapper, this, std::placeholders::_1, std::placeholders::_2),
        tsTime);

    std::cerr << "qow final = " << qow(79,67) << std::endl;
}

void IceGrowth::initializeThicknesses()
{
    cice = cice0;
    overElements(std::bind(&IceGrowth::initializeThicknessesElement, this, std::placeholders::_1,
                     std::placeholders::_2),
        TimestepTime());
}

// Divide by ice concentration to go from cell-averaged to ice-averaged values,
// but only if ice concentration is non-zero.
void IceGrowth::initializeThicknessesElement(size_t i, const TimestepTime&)
{
    deltaCIce[i] = 0;

    if (cice0[i] > 0 && hIceCell[i] > 0) {
        hice[i] = hice0[i] = hIceCell[i] / cice0[i];
        hsnow[i] = hsnow0[i] = hSnowCell[i] / cice0[i];
    } else {
        hice[i] = hice0[i] = 0.;
        hsnow[i] = hsnow0[i] = 0.;
        cice[i] = 0.;
    }
}

void IceGrowth::newIceFormation(size_t i, const TimestepTime& tst)
{
    bool doPrint = (i == ModelArray::indexFromLocation(ModelArray::Type::H, {79,67}));
    if (doPrint) std::cerr << "before new ice: qow = " << qow[i] << " hice = " << hice[i] << std::endl;
    // Flux cooling the ocean from open water
    // TODO Add assimilation fluxes here
    double coolingFlux = qow[i];
    // Temperature change of the mixed layer during this timestep
    double deltaTml = -coolingFlux / mixedLayerBulkHeatCapacity[i] * tst.step;
    // Initial temperature
    double t0 = sst[i];
    // Freezing point temperature
    double tf0 = tf[i];
    // Final temperature
    double t1 = t0 + deltaTml;

    if (doPrint) std::cerr << "cooling flux = " << coolingFlux << " T1 = " << t1 << " Tf = " << tf0;
    // deal with cooling below the freezing point
    if (t1 < tf0) {
        // Heat lost cooling the mixed layer to freezing point
        double sensibleFlux = (tf0 - t0) / deltaTml * coolingFlux;
        // Any heat beyond that is latent heat forming new ice
        double latentFlux = coolingFlux - sensibleFlux;

        qow[i] = sensibleFlux;
        newice[i] = latentFlux * tst.step * (1 - cice[i]) / (Ice::Lf * Ice::rho);

    if (doPrint) std::cerr << " qs = " << sensibleFlux << " ql = " << latentFlux << " qow = " << qow[i] << " newice = " << newice[i];
    } else {
        newice[i] = 0;
    }
    if (doPrint) std::cerr << std::endl << "after new ice: qow = " << qow[i] << " hice = " << hice[i] << std::endl;
}

// Update thickness with concentration
static double updateThickness(double& thick, double newConc, double deltaC, double deltaV)
{
    return thick += (deltaV - thick * deltaC) / newConc;
}

void IceGrowth::lateralIceSpread(size_t i, const TimestepTime& tstep)
{
    bool doPrint = (i == ModelArray::indexFromLocation(ModelArray::Type::H, {79,67}));
    if (doPrint) std::cerr << "before ice spread: qow = " << qow[i] << " hice = " << hice[i] << std::endl;
    deltaCMelt[i] = 0;
    deltaCFreeze[i] = 0;
    iLateral->freeze(
        tstep, hice[i], hsnow[i], deltaHi[i], newice[i], cice[i], qow[i], deltaCFreeze[i]);
    if (deltaHi[i] < 0) {
        // Note that the cell-averaged hice0 is converted to a ice averaged value
        iLateral->melt(tstep, hice0[i], hsnow[i], deltaHi[i], cice[i], qow[i], deltaCMelt[i]);
    }
    deltaCIce[i] = deltaCFreeze[i] + deltaCMelt[i];
    cice[i] = (hice[i] > 0 || newice[i] > 0) ? cice[i] + deltaCIce[i] : 0;
    if (cice[i] >= IceMinima::c()) {
        // The updated ice thickness must conserve volume
        updateThickness(hice[i], cice[i], deltaCIce[i], newice[i]);
        if (deltaCIce[i] < 0) {
            // Snow is lost if the concentration decreases, and energy is returned to the ocean
            qow[i] -= deltaCIce[i] * hsnow[i] * Water::Lf * Ice::rhoSnow / tstep.step;
        } else {
            // Update snow thickness. Currently no new snow is implemented
            updateThickness(hsnow[i], cice[i], deltaCIce[i], 0);
        }
    }
    if (doPrint) std::cerr << "after ice spread: qow = " << qow[i] << " hice = " << hice[i] << std::endl;

}

void IceGrowth::applyLimits(size_t i, const TimestepTime& tstep)
{
    bool doPrint = (i == ModelArray::indexFromLocation(ModelArray::Type::H, {79,67}));
    if (doPrint) std::cerr << "cice = " << cice[i] << " hice = " << hice[i] << std::endl;
    if ((0. < cice[i] && cice[i] < IceMinima::c()) || (0. < hice[i] && hice[i] < IceMinima::h())) {
        qow[i] += cice[i] * Water::Lf * (hice[i] * Ice::rho + hsnow[i] * Ice::rhoSnow) / tstep.step;
        hice[i] = 0;
        cice[i] = 0;
        hsnow[i] = 0;
    }
    if (hice[i] > 6.371e6)
        std::cerr << "Suspected polychoronya at " << i << std::endl;
    if (doPrint) std::cerr << "qow after apply limits = " << qow[i] << std::endl;

}
} /* namespace Nextsim */
