/*!
 * @file PrognosticData.cpp
 * @date Sep 8, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/PrognosticData.hpp"
#include "include/IFreezingPoint.hpp"

namespace Nextsim {
double PrognosticData::m_dt;
IFreezingPoint PrognosticData::m_freezer;
} /* namespace Nextsim */
