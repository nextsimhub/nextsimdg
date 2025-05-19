/*!
 * @file PrognosticData.cpp
 *
 * @date 21 Nov 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#include "include/PrognosticData.hpp"

#include "include/Finalizer.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/NextsimModule.hpp"
#include "include/gridNames.hpp"

namespace Nextsim {

PrognosticData::PrognosticData()
    : m_dt(1)
    , m_snow(ModelArray::Type::H)
    , m_tice(ModelArray::Type::Z)
    , m_damage(ModelArray::Type::H)
    , hiceAdvection(ModelArray::AdvectionType)
    , ciceAdvection(ModelArray::AdvectionType)
    , pAtmBdy(0)
    , pOcnBdy(0)
    , pDynamics(0)

{
    getStore().registerArray(Protected::H_ICE, &hiceAdvection, RO);
    getStore().registerArray(Protected::C_ICE, &ciceAdvection, RO);
    getStore().registerArray(Protected::H_SNOW, &m_snow, RO);
    getStore().registerArray(Protected::T_ICE, &m_tice, RO);
    getStore().registerArray(Protected::DAMAGE, &m_damage, RO);
    getStore().registerArray(Shared::H_ICE_DG, &hiceAdvection, RW);
    getStore().registerArray(Shared::C_ICE_DG, &ciceAdvection, RW);
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

void copyAllComponents(const ModelArray& source, ModelArray& sink)
{
    if (source.nComponents() == sink.nComponents()) {
        sink = source;
    } else if (source.nComponents() == 1) {
        sink.component(0) = source.data();
    } else {
        std::string err = std::string("PrognosticData::copyAllComponents: Expected 1 or ")
            + std::to_string(sink.nComponents()) + " components, got "
            + std::to_string(source.nComponents()) + " components.";
        throw std::runtime_error(err);
    }
}

void PrognosticData::setData(const ModelState::DataMap& ms)
{

    if (ms.count(maskName)) {
        setOceanMask(ms.at(maskName));
    } else {
        noLandMask();
    }

    copyMeanComponent(ms.at(ticeName), m_tice);
    copyMeanComponent(ms.at(hsnowName), m_snow);
    // Damage is an optional field, and defaults to 1, if absent
    if (ms.count(damageName) > 0) {
        copyMeanComponent(ms.at(damageName), m_damage);
    } else {
        m_damage.resize();
        m_damage = 1.;
    }

    // Copy the full DG data
    hiceAdvection = 0;
    ciceAdvection = 0;
    copyAllComponents(ms.at(hiceName), hiceAdvection);
    copyAllComponents(ms.at(ciceName), ciceAdvection);

    pAtmBdy->setData(ms);
    pOcnBdy->setData(ms);
    pDynamics->setData(ms);
    iceGrowth.setData(ms);
}

void PrognosticData::update(const TimestepTime& tst)
{
    pOcnBdy->updateBefore(tst);
    pAtmBdy->update(tst);

    // Take the updated values of the true ice and snow thicknesses, and reset hice0 and hsnow0
    // IceGrowth updates its own fields during update
    iceGrowth.update(tst);
    updatePrognosticFields();

    pDynamics->update(tst);

    updateDynamicsFields();

    pOcnBdy->updateAfter(tst);
}

void PrognosticData::updatePrognosticFields()
{
    ModelArrayRef<Shared::H_ICE, RO> hiceTrueUpd(getStore());
    ModelArrayRef<Shared::C_ICE, RO> ciceUpd(getStore());
    ModelArrayRef<Shared::H_SNOW, RO> hsnowTrueUpd(getStore());
    ModelArrayRef<Shared::T_ICE, RO> ticeUpd(getStore());
    ModelArrayRef<Shared::DAMAGE, RO> damageUpd(getStore());

    // Calculate the cell average thicknesses
    HField hiceUpd = hiceTrueUpd * ciceUpd;
    HField hsnowUpd = hsnowTrueUpd * ciceUpd;

    // Update the DG0 component of the DG fields
    hiceAdvection.component(0) = hiceUpd.data();
    ciceAdvection.component(0) = ciceUpd.allComponents();
    m_snow.setData(hsnowUpd);
    m_tice.setData(ticeUpd);
    m_damage.setData(damageUpd);
}

void PrognosticData::updateDynamicsFields()
{
    ModelArrayRef<Shared::H_SNOW, RO> hsnowTrueUpd(getStore());
    ModelArrayRef<Shared::T_ICE, RO> ticeUpd(getStore());
    ModelArrayRef<Shared::DAMAGE, RO> damageUpd(getStore());

    // Calculate the cell average thicknesses
    HField hsnowUpd;
    hsnowUpd.setData(hsnowTrueUpd.allComponents() * ciceAdvection.component(0));

    m_snow.setData(hsnowUpd);
    m_tice.setData(ticeUpd);
    m_damage.setData(damageUpd);
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
                             { "mask", ModelArray(oceanMask()) }, // make a copy
                             { "hice", hiceAdvection },
                             { "cice", ciceAdvection },
                             { "hsnow", mask(m_snow) },
                             { "tice", mask(m_tice) },
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
