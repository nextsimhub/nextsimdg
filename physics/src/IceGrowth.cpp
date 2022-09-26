/*!
 * @file IceGrowth.cpp
 *
 * @date Mar 15, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/IceGrowth.hpp"

#include "include/constants.hpp"

namespace Nextsim {

double IceGrowth::minc;
double IceGrowth::minh;

static const double mincDefault = 1e-12;
static const double minhDefault = 0.01;

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
        { keyMap.at(MINC_KEY), ConfigType::NUMERIC, { "0", "1" }, std::to_string(mincDefault), "",
            "Minimum allowed ice concentration." },
        { keyMap.at(MINH_KEY), ConfigType::NUMERIC, { "0", "∞" }, std::to_string(minhDefault), "m",
            "Minimum allowed ice thickness." },
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
    minc = Configured::getConfiguration(keyMap.at(MINC_KEY), mincDefault);
    minh = Configured::getConfiguration(keyMap.at(MINH_KEY), minhDefault);

    // Configure the vertical and lateral growth modules
    iVertical = std::move(Module::getInstance<IIceThermodynamics>());
    iLateral = std::move(Module::getInstance<ILateralIceSpread>());
    tryConfigure(*iVertical);
    tryConfigure(*iLateral);
}

ConfigMap IceGrowth::getConfiguration() const
{
    return {
        { keyMap.at(MINC_KEY), minc },
        { keyMap.at(MINH_KEY), minh },
    };
}

void IceGrowth::update(const TimestepTime& tsTime)
{

    // Copy the ice data from the prognostic fields to the modifiable fields.
    cice = cice0;
    overElements(std::bind(&IceGrowth::initializeThicknesses, this, std::placeholders::_1,
                     std::placeholders::_2),
        tsTime);

    iVertical->update(tsTime);
    // new ice formation
    overElements(
        std::bind(&IceGrowth::updateWrapper, this, std::placeholders::_1, std::placeholders::_2),
        tsTime);
}

// Divide by ice concentration to go from cell-averaged to ice-averaged values,
// but only if ice concentration is non-zero.
void IceGrowth::initializeThicknesses(size_t i, const TimestepTime&)
{
    if (cice0[i] > 0) {
        hice[i] = hice0[i] = hIceCell[i] / cice0[i];
        hsnow[i] = hsnow0[i] = hSnowCell[i] / cice0[i];
    } else {
        hice[i] = hice0[i] = 0.;
        hsnow[i] = hsnow0[i] = 0.;
    }
}

void IceGrowth::newIceFormation(size_t i, const TimestepTime& tst)
{
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

    // deal with cooling below the freezing point
    if (t1 < tf0) {
        // Heat lost cooling the mixed layer to freezing point
        double sensibleFlux = (tf0 - t0) / deltaTml * coolingFlux;
        // Any heat beyond that is latent heat forming new ice
        double latentFlux = coolingFlux - sensibleFlux;

        qow[i] = sensibleFlux;
        newice[i] = latentFlux * tst.step * (1 - cice[i]) / (Ice::Lf * Ice::rho);
    } else {
        newice[0] = 0;
    }
}

// Update thickness with concentration
static double updateThickness(double& thick, double newConc, double deltaC, double deltaV)
{
    return thick += (deltaV - thick * deltaC) / newConc;
}

void IceGrowth::lateralIceSpread(size_t i, const TimestepTime& tstep)
{
    iLateral->freeze(
        tstep, hice[i], hsnow[i], deltaHi[i], newice[i], cice[i], qow[i], deltaCFreeze[i]);
    if (deltaHi[i] < 0) {
        // Note that the cell-averaged hice0 is converted to a ice averaged value
        iLateral->melt(tstep, hice0[i], hsnow[i], deltaHi[i], cice[i], qow[i], deltaCMelt[i]);
    }
    double deltaC = deltaCFreeze[i] + deltaCMelt[i];
    cice[i] += deltaC;
    if (cice[i] >= minc) {
        // The updated ice thickness must conserve volume
        updateThickness(hice[i], cice[i], deltaC, newice[i]);
        if (deltaC < 0) {
            // Snow is lost if the concentration decreases, and energy is returned to the ocean
            qow[i] -= deltaC * hsnow[i] * Water::Lf * Ice::rhoSnow / tstep.step;
        } else {
            // Update snow thickness. Currently no new snow is implemented
            updateThickness(hsnow[i], cice[i], deltaC, 0);
        }
    }
}

void IceGrowth::applyLimits(size_t i, const TimestepTime& tstep)
{
    if (cice[i] < minc || hice[i] < minh) {
        qow[i] += cice[i] * Water::Lf * (hice[i] * Ice::rho + hsnow[i] * Ice::rhoSnow) / tstep.step;
        hice[i] = 0;
        cice[i] = 0;
        hsnow[i] = 0;
    }
}
} /* namespace Nextsim */
