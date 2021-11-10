/*!
 * @file PrognosticData.cpp
 * @date Sep 8, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/PrognosticData.hpp"
#include "include/IFreezingPoint.hpp"
#include "ModuleLoader.hpp"
namespace Nextsim {

double PrognosticData::m_dt = 0;
IFreezingPoint* PrognosticData::m_freezer = nullptr;

PrognosticData::PrognosticData()
    : m_conc(0)
    , m_snow(0)
    , m_sss(0)
    , m_sst(0)
    , m_thick(0)
{
    if (!m_freezer) m_freezer = ModuleLoader::getLoader().getImplementation<IFreezingPoint>();
}
} /* namespace Nextsim */
