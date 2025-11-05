/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Einar Ólason <einar.olason@nersc.no>
 */

#include "include/PrognosticData.hpp"

#include "include/FieldAdvection.hpp"
#include "include/Finalizer.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/NextsimModule.hpp"
#include "include/gridNames.hpp"

namespace Nextsim {

static const std::string checkFieldsKey = "debug.check_fields";
static const std::string checkFieldsFastKey = "debug.check_fields_fast";
static const std::map<int, std::string> keyMap
    = { { PrognosticData::CHECKFIELDS_KEY, checkFieldsKey },
          { PrognosticData::CHECKFIELDSFAST_KEY, checkFieldsFastKey } };
static constexpr bool checkFieldsDefault = false;
static constexpr bool checkFieldsFastDefault = true;

PrognosticData::PrognosticData()
    : m_dt(1)
    , hice(ModelArray::AdvectionType, { 0, 50 })
    , cice(ModelArray::AdvectionType, { 0, 1 })
    , damage(ModelArray::AdvectionType, { 0, 1 })
    , hsnow(ModelArray::AdvectionType, { 0, 10 })
    , pAtmBdy(nullptr)
    , pOcnBdy(nullptr)
    , pDynamics(nullptr)
{
    getStore().registerArray(Shared::DAMAGE, &damage, RW);
    getStore().registerArray(Shared::H_ICE_DG, &hice, RW);
    getStore().registerArray(Shared::C_ICE_DG, &cice, RW);
    getStore().registerArray(Shared::H_SNOW_DG, &hsnow, RW);
}

void PrognosticData::configure()
{
    // Register finalizers before calling configure.
    Finalizer::registerUnique(Module::finalize<IAtmosphereBoundary>);
    Finalizer::registerUnique(Module::finalize<IOceanBoundary>);
    Finalizer::registerUnique(Module::finalize<IDynamics>);

    pAtmBdy = &Module::getImplementation<IAtmosphereBoundary>();
    tryConfigure(pAtmBdy);

    pOcnBdy = &Module::getImplementation<IOceanBoundary>();
    tryConfigure(pOcnBdy);

    pDynamics = &Module::getImplementation<IDynamics>();
    tryConfigure(pDynamics);

    tryConfigure(iceGrowth);

    checkAll() = Configured::getConfiguration(keyMap.at(CHECKFIELDS_KEY), checkFieldsDefault);
    checkFast
        = Configured::getConfiguration(keyMap.at(CHECKFIELDSFAST_KEY), checkFieldsFastDefault);
    if (checkAll()) {
        for (const auto& field : getStore().getAllData()) {
            addChecks({ { field.first, field.second } });
        }
    } else if (checkFast) {
        addChecks({
            { "thickness", &hice },
            { "concentration", &cice },
        });
    }
}

// Copies an HField from a source ModelArray that is either an HField or a DGField.
void copyMeanComponent(const ModelArray& source, ModelArray& sink)
{
    if (source.nComponents() > 1) {
        sink.setData(source.data().col(0));
    } else {
        sink = source;
    }
}

void PrognosticData::setData(const ModelState::DataMap& ms)
{
    if (ms.count(maskName)) {
        setOceanMask(ms.at(maskName));
    } else {
        noLandMask();
    }

    // Copy the full DG data
    hice = 0;
    cice = 0;
    hsnow = 0;
    damage = 0;
    hice = ms.at(hiceName);
    cice = ms.at(ciceName);
    hsnow = ms.at(hsnowName);
    // Damage is an optional field, and defaults to 1 in the mean field, 0 in higher components
    // if absent.
    if (ms.count(damageName) > 0) {
        damage = ms.at(damageName);
    } else {
        damage.component(0) = 1.;
    }

    pAtmBdy->setData(ms);
    pOcnBdy->setData(ms);
    pDynamics->setData(ms);
    iceGrowth.setData(ms);
}

void PrognosticData::update(const TimestepTime& tst)
{
    // Prepare everything
    pOcnBdy->updateBefore(tst);
    pAtmBdy->update(tst);
    pDynamics->prepareAdvection();

    // Take the updated values of the true ice and snow thicknesses, and reset hice0 and hsnow0
    // IceGrowth updates its own fields during update
    iceGrowth.update(tst);

    // Dynamics
    pDynamics->update(tst);

    // Update the ocean after ice growth (or send fields to the coupler)
    pOcnBdy->updateAfter(tst);

    try {
        checkFields();
    } catch (const std::exception& e) {
        throw std::runtime_error("PrognosticData::update: " + std::string(e.what()));
    }
}

// Gets the diagnostic data from all subcomponents
ModelState PrognosticData::getStateDiagnostic() const
{
    ModelState state = getStatePrognostic();

    // Get the prognostic data from the dynamics, including the full dynamics state
    state.merge(pDynamics->getStateDiagnostic());
    state.merge(iceGrowth.getStateDiagnostic());
    state.merge(pAtmBdy->getStateDiagnostic());
    state.merge(pOcnBdy->getStateDiagnostic());

    return state;
}

// Gets the prognostic data from all subcomponents
ModelState PrognosticData::getStatePrognostic() const
{
    ModelState state = { {
                             { maskName, ModelArray(oceanMask()) }, // make a copy
                             { hiceName, hice },
                             { ciceName, cice },
                             { hsnowName, hsnow },
                         },
        ModelComponent::getConfiguration() };

    // Get the prognostic data from the dynamics, including the full dynamics state
    state.merge(pDynamics->getStatePrognostic());
    state.merge(iceGrowth.getStatePrognostic());
    state.merge(pAtmBdy->getStatePrognostic());
    state.merge(pOcnBdy->getStatePrognostic());

    return state;
}

PrognosticData::HelpMap& PrognosticData::getHelpText(HelpMap& map, bool getAll)
{
    map["debug"] = {
        { checkFieldsKey, ConfigType::BOOLEAN, { "true", "false" },
            ConfigurationHelp::toString(checkFieldsDefault), "",
            "Set to true to check if all variables in the ModelArrayReferenceStore fall within a "
            "reasonable physical range." },
        { checkFieldsFastKey, ConfigType::BOOLEAN, { "true", "false" },
            ConfigurationHelp::toString(checkFieldsFastDefault), "",
            "Set to true to check if thickness, concentration, and velocities fall within a "
            "reasonable physical range." }
    };

    return map;
}
PrognosticData::HelpMap& PrognosticData::getHelpRecursive(HelpMap& map, bool getAll)
{
    Module::getHelpRecursive<IAtmosphereBoundary>(map, getAll);
    Module::getHelpRecursive<IOceanBoundary>(map, getAll);
    Module::getHelpRecursive<IDynamics>(map, getAll);
    IceGrowth::getHelpRecursive(map, getAll);
    getHelpText(map, getAll);
    return map;
}

} /* namespace Nextsim */
