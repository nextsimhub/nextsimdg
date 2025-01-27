/*!
 * @file PrognosticData.cpp
 *
 * @date 27 Jan 2025
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
    , m_thick(ModelArray::Type::H)
    , m_conc(ModelArray::Type::H)
    , m_snow(ModelArray::Type::H)
    , m_tice(ModelArray::Type::Z)
    , m_damage(ModelArray::Type::H)
    , pAtmBdy(0)
    , pOcnBdy(0)
    , pDynamics(0)

{
    getStore().registerArray(Shared::C_ICE, &m_conc, RW);
    getStore().registerArray(Shared::DAMAGE, &m_damage, RW);
    getStore().registerArray(Shared::H_ICE, &m_thick, RW);
    getStore().registerArray(Shared::H_SNOW, &m_snow, RW);
    getStore().registerArray(Shared::T_ICE, &m_tice, RW);
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

    copyMeanComponent(ms.at(hiceName), m_thick);
    copyMeanComponent(ms.at(ciceName), m_conc);
    copyMeanComponent(ms.at(ticeName), m_tice);
    copyMeanComponent(ms.at(hsnowName), m_snow);
    // Damage is an optional field, and defaults to 1, if absent
    if (ms.count(damageName) > 0) {
        copyMeanComponent(ms.at(damageName), m_damage);
    } else {
        m_damage.resize();
        m_damage = 1.;
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

    // Take the updated values of the true ice and snow thicknesses, and reset hice0 and hsnow0
    // IceGrowth updates its own fields during update
    iceGrowth.update(tst);

    pDynamics->update(tst);

    pOcnBdy->updateAfter(tst);
}

// Gets all of the prognostic data, including that in the dynamics
ModelState PrognosticData::getState() const
{
    ModelArrayRef<Shared::SSS> sss(getStore());
    ModelArrayRef<Shared::SST> sst(getStore());

    // Get the prognostic data from the dynamics, including the full dynamics state
    ModelState dynamicsState = pDynamics->getState();
    // clang-format off
    ModelState localState = { {
                 { "mask", ModelArray(oceanMask()) }, // make a copy
                 { "hice", mask(m_thick) },
                 { "cice", mask(m_conc) },
                 { "hsnow", mask(m_snow) },
                 { "tice", mask(m_tice) },
                 { "sst", mask(sst.data()) },
                 { "sss", mask(sss.data()) },
             },
        {} };
    // clang-format on

    // Use the dynamics values of any duplicated fields
    ModelState state(dynamicsState);
    state.merge(localState);

    return state;
}

// Recursively gets the data from all subcomponents
ModelState PrognosticData::getStateRecursive(const OutputSpec& os) const
{
    ModelState state;
    /* If allComponents is set on the OutputSpec, then for any duplicate fields, the subsystems
     * take priority, otherwise the fields held by PrognosticData itself. Note that std::map::merge
     * will not overwrite existing keys, so the first one that exists will survive.
     */
    if (os.allComponents()) {
        state.merge(pAtmBdy->getStateRecursive(os));
        state.merge(iceGrowth.getStateRecursive(os));
        state.merge(pDynamics->getStateRecursive(os));
        state.merge(getState());
    } else {
        state.merge(getState());
        state.merge(pAtmBdy->getStateRecursive(os));
        state.merge(iceGrowth.getStateRecursive(os));
        state.merge(pDynamics->getStateRecursive(os));
    }
    // OceanBoundary does not contribute to the output model state
    return os ? state : ModelState();
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
