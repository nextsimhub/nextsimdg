/*!
 * @file PrognosticData.cpp
 *
 * @date Mar 1, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/PrognosticData.hpp"

#include "include/ModelArrayRef3.hpp"
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

{
    getStore().registerArray(Protected::H_ICE, &m_thick);
    getStore().registerArray(Protected::C_ICE, &m_conc);
    getStore().registerArray(Protected::H_SNOW, &m_snow);
    getStore().registerArray(Protected::T_ICE, &m_tice);
}

void PrognosticData::configure()
{
    tryConfigure(iceGrowth);

    pAtmBdy = &Module::getImplementation<IAtmosphereBoundary>();
    tryConfigure(pAtmBdy);

    pOcnBdy = &Module::getImplementation<IOceanBoundary>();
    tryConfigure(pOcnBdy);
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
    if (ms.count("u")) {
        m_u = ms.at("u");
    }
    if (ms.count("v")) {
        m_v = ms.at("v");
    }

    iceGrowth.setData(ms);
    pAtmBdy->setData(ms);
    pOcnBdy->setData(ms);
}

void PrognosticData::update(const TimestepTime& tst)
{
    pAtmBdy->update(tst);
    pOcnBdy->updateBefore(tst);
    iceGrowth.update(tst);
    pOcnBdy->updateAfter(tst);

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
    return { {
                 { "mask", ModelArray(oceanMask()) }, // make a copy
                 { "hice", mask(m_thick) },
                 { "cice", mask(m_conc) },
                 { "hsnow", mask(m_snow) },
                 { "tice", mask(m_tice) },
             },
        {} };
}

ModelState PrognosticData::getStateRecursive(const OutputSpec& os) const
{
    ModelState state(getState());
    state.merge(iceGrowth.getStateRecursive(os));
    state.merge(pAtmBdy->getStateRecursive(os));
    return os ? state : ModelState();
}

PrognosticData::HelpMap& PrognosticData::getHelpText(HelpMap& map, bool getAll) { return map; }
PrognosticData::HelpMap& PrognosticData::getHelpRecursive(HelpMap& map, bool getAll)
{
    IceGrowth::getHelpRecursive(map, getAll);
    Module::getHelpRecursive<IAtmosphereBoundary>(map, getAll);
    Module::getHelpRecursive<IOceanBoundary>(map, getAll);
    return map;
}

} /* namespace Nextsim */
