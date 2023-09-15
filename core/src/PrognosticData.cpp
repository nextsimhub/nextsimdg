/*!
 * @file PrognosticData.cpp
 *
 * @date 7 Sep 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/PrognosticData.hpp"

#include "include/ModelArrayRef.hpp"
#include "include/Module.hpp"

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
    getStore().registerArray(Protected::H_ICE, &m_thick);
    getStore().registerArray(Protected::C_ICE, &m_conc);
    getStore().registerArray(Protected::H_SNOW, &m_snow);
    getStore().registerArray(Protected::T_ICE, &m_tice);
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

    // Calculate the cell average thicknesses
    HField hiceUpd = hiceTrueUpd * ciceUpd;
    HField hsnowUpd = hsnowTrueUpd * ciceUpd;

    m_thick.setData(hiceUpd);
    m_conc.setData(ciceUpd);
    m_snow.setData(hsnowUpd);
    m_tice.setData(ticeUpd);
}

ModelState PrognosticData::getState() const
{
    ModelArrayRef<Protected::SST> sst(getStore());
    ModelArrayRef<Protected::SSS> sss(getStore());
    return { {
                 { "mask", ModelArray(oceanMask()) }, // make a copy
                 { "hice", mask(m_thick) },
                 { "cice", mask(m_conc) },
                 { "hsnow", mask(m_snow) },
                 { "tice", mask(m_tice) },
                 { "sst", mask(sst.data()) },
                 { "sss", mask(sss.data()) },
             },
        {} };
}

ModelState PrognosticData::getStateRecursive(const OutputSpec& os) const
{
    ModelState state(getState());
    state.merge(pAtmBdy->getStateRecursive(os));
    state.merge(iceGrowth.getStateRecursive(os));
    // Neither OceanBdoudary nor Dynamics contribute to the output model state
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
