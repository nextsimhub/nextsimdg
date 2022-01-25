/*!
 * @file PrognosticData.cpp
 * @date Sep 8, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/PrognosticData.hpp"
#include "include/IFreezingPoint.hpp"
#include "include/ModuleLoader.hpp"
namespace Nextsim {

double PrognosticData::m_dt = 0;
IFreezingPoint* PrognosticData::m_freezer = nullptr;

PrognosticData::PrognosticData()
    : PrognosticData(1)
{}

PrognosticData::PrognosticData(int nIceLayers)
    : m_conc(0)
    , m_snow(0)
    , m_sss(0)
    , m_sst(0)
    , m_thick(0)
    , m_tice(nIceLayers, 0.)
{
//    m_tice.resize(nIceLayers);
}

void PrognosticData::configure()
{
    ModuleLoader& loader = ModuleLoader::getLoader();
    m_freezer = &loader.getImplementation<IFreezingPoint>();
    tryConfigure(m_freezer);
}
} /* namespace Nextsim */
