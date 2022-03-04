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
    registerSharedArray(SharedArray::H_ICE, &m_thick);
    registerSharedArray(SharedArray::C_ICE, &m_conc);
    registerSharedArray(SharedArray::H_SNOW, &m_snow);
    registerSharedArray(SharedArray::T_ICE, &m_tice);
}

PrognosticData::~PrognosticData()
{
}

} /* namespace Nextsim */
