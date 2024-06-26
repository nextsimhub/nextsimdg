/*!
 * @file PrognosticData.cpp
 *
 * @date 7 Sep 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/PrognosticData.hpp"

#include "include/ModelArrayRef.hpp"
#include "include/Module.hpp"
#include "include/gridNames.hpp"

namespace Nextsim {

PrognosticData::PrognosticData()
    : m_dt(1)
    , m_thick(ModelArray::Type::H)
    , m_conc(ModelArray::Type::H)
    , m_snow(ModelArray::Type::H)
    , m_tice(ModelArray::Type::Z)
    , pAtmBdy(0)
    , pOcnBdy(0)
    , pDynamics(0)

{
    getStore().registerArray(Protected::H_ICE, &m_thick, RO);
    getStore().registerArray(Protected::C_ICE, &m_conc, RO);
    getStore().registerArray(Protected::H_SNOW, &m_snow, RO);
    getStore().registerArray(Protected::T_ICE, &m_tice, RO);
    getStore().registerArray(Protected::DAMAGE, &m_damage, RO);
}

void PrognosticData::configure()
{
    pAtmBdy = &Module::getImplementation<IAtmosphereBoundary>();
    tryConfigure(pAtmBdy);

    pOcnBdy = &Module::getImplementation<IOceanBoundary>();
    tryConfigure(pOcnBdy);

    pDynamics = &Module::getImplementation<IDynamics>();
    tryConfigure(pDynamics);

    tryConfigure(iceGrowth);
}

void PrognosticData::setData(const ModelState::DataMap& ms)
{

    if (ms.count("mask")) {
        setOceanMask(ms.at("mask"));
    } else {
        noLandMask();
    }

    m_thick = ms.at("hice");
    m_conc = ms.at("cice");
    m_tice = ms.at("tice");
    m_snow = ms.at("hsnow");
    // Damage is an optional field, and defaults to zero, if absent
    if (ms.count(damageName) > 0) {
        m_damage = ms.at(damageName);
    } else {
        m_damage.resize();
        m_damage = 0.5;
    }

    pAtmBdy->setData(ms);
    pOcnBdy->setData(ms);
    pDynamics->setData(ms);
    iceGrowth.setData(ms);
}

void PrognosticData::update(const TimestepTime& tst)
{
    ModelArrayRef<Shared::T_ICE, RW> ticeUpd(getStore());

    pOcnBdy->updateBefore(tst);
    pAtmBdy->update(tst);

    // Fill the values of the true ice and snow thicknesses.
    iceGrowth.initializeThicknesses();
    // Fill the updated ice temperature array
    ticeUpd.data().setData(m_tice);
    pDynamics->update(tst);
    updatePrognosticFields();

    // Take the updated values of the true ice and snow thicknesses, and reset hice0 and hsnow0
    // IceGrowth updates its own fields during update
    iceGrowth.update(tst);
    updatePrognosticFields();

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

    m_thick.setData(hiceUpd);
    m_conc.setData(ciceUpd);
    m_snow.setData(hsnowUpd);
    m_tice.setData(ticeUpd);
    m_damage.setData(damageUpd);
}

ModelState PrognosticData::getState() const
{
    ModelArrayRef<Protected::SST> sst(getStore());
    ModelArrayRef<Protected::SSS> sss(getStore());
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
    // Get the state from the dynamics (ice velocity). This allows the
    // dynamics to define its own dimensions for the velocity grid.
    localState.merge(pDynamics->getState());

    // Merge in the damage field, if the dynamics uses it.
    if (pDynamics->usesDamage()) {
        ModelState damageState = { {
                                       { "damage", mask(m_damage) },
                                   },
            {} };
        localState.merge(damageState);
    }
    return localState;
}

ModelState PrognosticData::getStateRecursive(const OutputSpec& os) const
{
    ModelState state;
    /* If allComponents is set on the OutputSpec, then for any duplicate fields, the subsystems
     * take priority, otherwise the fields held by PrognosticData itself. Note that merge will not
     * overwrite existing keys, so the first one that exists will survive.
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
