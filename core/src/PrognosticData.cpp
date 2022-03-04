/*!
 * @file PrognosticData.cpp
 *
 * @date Mar 1, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/PrognosticData.hpp"

namespace Nextsim {

PrognosticData::PrognosticData()
    : m_dt(1)
{
    registerProtectedArray(ProtectedArray::H_ICE, &m_thick);
    registerProtectedArray(ProtectedArray::C_ICE, &m_conc);
    registerProtectedArray(ProtectedArray::H_SNOW, &m_snow);
    registerProtectedArray(ProtectedArray::T_ICE, &m_tice);
}

void PrognosticData::setData(const ModelState& ms)
{
    m_thick = *ms["hice"];
    m_conc = *ms["cice"];
    m_tice = *ms["tice"];
    m_snow = *ms["hsnow"];
    m_u = *ms["u"];
    m_v = *ms["v"];
}

ModelState PrognosticData::getState() const
{
    ModelState ms;

    return ms;
} /* namespace Nextsim */
