/*!
 * @file PrognosticData.cpp
 * @date Sep 8, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/PrognosticData.hpp"
#include "include/IFreezingPoint.hpp"
#include "include/Module.hpp"
namespace Nextsim {

double PrognosticData::m_dt = 0;
IFreezingPoint* PrognosticData::m_freezer = nullptr;

PrognosticData::PrognosticData()
    : PrognosticData(1)
{
}

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

PrognosticData::PrognosticData(const PrognosticGenerator& up)
    : m_conc(up.updatedIceConcentration())
    , m_snow(up.updatedSnowThickness())
    , m_thick(up.updatedIceThickness())
    , m_tice()
    , m_sst(up.seaSurfaceTemperature())
    , m_sss(up.seaSurfaceSalinity())
{
    m_tice = up.updatedIceTemperatures();

    configure();
}

PrognosticData& PrognosticData::operator=(const PrognosticGenerator& up)
{
    m_thick = up.updatedIceThickness();
    m_conc = up.updatedIceConcentration();
    m_snow = up.updatedSnowThickness();

    m_tice = up.updatedIceTemperatures();

    m_sst = up.seaSurfaceTemperature();
    m_sss = up.seaSurfaceSalinity();

    configure();

    return *this;
}

void PrognosticData::configure()
{
    m_freezer = &Module::getImplementation<IFreezingPoint>();
    tryConfigure(m_freezer);
}

PrognosticData& PrognosticData::updateAndIntegrate(const IPrognosticUpdater& updater)
{
    m_thick = updater.updatedIceThickness();
    m_conc = updater.updatedIceConcentration();
    m_snow = updater.updatedSnowThickness();

    copyInIceLayerData(updater.updatedIceTemperatures(), m_tice);
    return *this;
}

PrognosticData& PrognosticData::setSeaSurface(double sst, double sss)
{
    m_sst = sst;
    m_sss = sss;
    return *this;
}

// Copy up to as many levels as there are currently.
// Fill missing layers with the lowest valid temperature
void PrognosticData::copyInIceLayerData(const std::vector<double>& src, std::vector<double>& tgt)
{
    int tLayers = tgt.size();
    int sLayers = src.size();

    tgt = src;
    if (tLayers != sLayers) {
        tgt.resize(tLayers);
        for (int i = sLayers; i < tLayers; ++i) {
            tgt[i] = tgt[sLayers - 1];
        }
    }
}
} /* namespace Nextsim */
