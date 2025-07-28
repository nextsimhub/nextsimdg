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

PrognosticData::PrognosticData()
    : m_dt(1)
    , hice(ModelArray::AdvectionType)
    , cice(ModelArray::AdvectionType)
    , damage(ModelArray::AdvectionType)
    , hsnow(ModelArray::AdvectionType)
    , pAtmBdy(0)
    , pOcnBdy(0)
    , pDynamics(0)

{
    getStore().registerArray(Protected::H_ICE, &hice, RO);
    getStore().registerArray(Protected::C_ICE, &cice, RO);
    getStore().registerArray(Protected::H_SNOW, &hsnow, RO);
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
    pOcnBdy->updateBefore(tst);
    pAtmBdy->update(tst);
    pDynamics->prepareAdvection();

    // Take the updated values of the true ice and snow thicknesses, and reset hice0 and hsnow0
    // IceGrowth updates its own fields during update
    iceGrowth.update(tst);
    updatePrognosticFields();

    pDynamics->update(tst);

    pOcnBdy->updateAfter(tst);
}

void PrognosticData::updatePrognosticFields()
{
    ModelArrayRef<Shared::H_ICE, RO> hiceTrueUpd(getStore());
    ModelArrayRef<Shared::C_ICE, RO> ciceUpd(getStore());
    ModelArrayRef<Shared::H_SNOW, RO> hsnowTrueUpd(getStore());
    ModelArrayRef<Shared::T_ICE, RO> ticeUpd(getStore());

    // Calculate the cell average thicknesses
    HField hiceUpd = hiceTrueUpd * ciceUpd;
    HField hsnowUpd = hsnowTrueUpd * ciceUpd;

    // Update the DG0 component of the DG fields
    hice.component(0) = hiceUpd.data();
    cice.component(0) = ciceUpd.allComponents();
    hsnow.component(0) = hsnowUpd.data();
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

PrognosticData::HelpMap& PrognosticData::getHelpText(HelpMap& map, bool getAll) { return map; }
PrognosticData::HelpMap& PrognosticData::getHelpRecursive(HelpMap& map, bool getAll)
{
    Module::getHelpRecursive<IAtmosphereBoundary>(map, getAll);
    Module::getHelpRecursive<IOceanBoundary>(map, getAll);
    Module::getHelpRecursive<IDynamics>(map, getAll);
    IceGrowth::getHelpRecursive(map, getAll);
    return map;
}

} /* namespace Nextsim */
