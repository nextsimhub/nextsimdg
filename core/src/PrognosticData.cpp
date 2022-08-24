/*!
 * @file PrognosticData.cpp
 *
 * @date Mar 1, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/PrognosticData.hpp"

#include "include/ModelArrayRef.hpp"

namespace Nextsim {

PrognosticData::PrognosticData()
    : m_dt(1)
    , m_thick(ModelArray::Type::H)
    , m_conc(ModelArray::Type::H)
    , m_snow(ModelArray::Type::H)
    , m_tice(ModelArray::Type::Z)

{
    registerProtectedArray(ProtectedArray::H_ICE, &m_thick);
    registerProtectedArray(ProtectedArray::C_ICE, &m_conc);
    registerProtectedArray(ProtectedArray::H_SNOW, &m_snow);
    registerProtectedArray(ProtectedArray::T_ICE, &m_tice);
}

void PrognosticData::configure() { tryConfigure(iceGrowth); }

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
}

void PrognosticData::update(const TimestepTime& tst)
{
    iceGrowth.update(tst);

    ModelArrayRef<SharedArray::H_ICE, RO> hiceTrueUpd;
    ModelArrayRef<SharedArray::C_ICE, RO> ciceUpd;
    ModelArrayRef<SharedArray::H_SNOW, RO> hsnowTrueUpd;
    ModelArrayRef<SharedArray::T_ICE, RO> ticeUpd;

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
    return {{
        { "mask", ModelArray(oceanMask()) }, // make a copy
        { "hice", mask(m_thick) },
        { "cice", mask(m_conc) },
        { "hsnow", mask(m_snow) },
        { "tice", mask(m_tice) },
    }, {}};
}

ModelState PrognosticData::getStateRecursive(const OutputSpec& os) const
{
    ModelState state(getState());
    state.merge(iceGrowth.getStateRecursive(os));
    return os ? state : ModelState();
}

PrognosticData::HelpMap& PrognosticData::getHelpText(HelpMap& map, bool getAll) { return map; }
PrognosticData::HelpMap& PrognosticData::getHelpRecursive(HelpMap& map, bool getAll)
{
    IceGrowth::getHelpRecursive(map, getAll);
    return map;
}

} /* namespace Nextsim */
