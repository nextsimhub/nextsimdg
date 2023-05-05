/*!
 * @file PrognosticData.cpp
 *
 * @date Mar 1, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/PrognosticData.hpp"

#include "include/ModelArrayRef.hpp"
#include "include/Module.hpp"

#include <iostream> // FIXME remove me
static const size_t itar = 64 * 128 + 64;

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
    registerProtectedArray(ProtectedArray::H_ICE, &m_thick);
    registerProtectedArray(ProtectedArray::C_ICE, &m_conc);
    registerProtectedArray(ProtectedArray::H_SNOW, &m_snow);
    registerProtectedArray(ProtectedArray::T_ICE, &m_tice);
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
    // Debugging MARs
    ModelArrayRef<ProtectedArray::HTRUE_ICE, MARBackingStore> hiceTrue0(getSharedArray());
    ModelArrayRef<SharedArray::H_ICE, MARBackingStore, RO> hiceTrueUpd(getSharedArray());
    ModelArrayRef<SharedArray::C_ICE, MARBackingStore, RO> ciceUpd(getSharedArray());
    // std::cerr << "PD::Update start: cice0=" << m_conc[itar] << " cice=" << ciceUpd[itar] << "
    // hice0=" << m_thick[itar] /*<< " htrue0=" << hiceTrue0[itar] */<< " htrue=" <<
    // hiceTrueUpd[itar] << std::endl;

    pOcnBdy->updateBefore(tst);
    pAtmBdy->update(tst);

    // The dynamics uses the SharedArrays of hice, cice and hsnow.
    // Fill the values of the true ice and snow thicknesses.
    // std::cerr << "Before initTh: cice0=" << m_conc[itar] << " cice=" << ciceUpd[itar] << "
    // hice0=" << m_thick[itar] /*<< " htrue0=" << hiceTrue0[itar] */<< " htrue=" <<
    // hiceTrueUpd[itar] << std::endl;
    iceGrowth.initializeThicknesses();
    // std::cerr << "After initTh: cice0=" << m_conc[itar] << " cice=" << ciceUpd[itar] << " hice0="
    // << m_thick[itar] /*<< " htrue0=" << hiceTrue0[itar] */<< " htrue=" << hiceTrueUpd[itar] <<
    // std::endl;
    pDynamics->update(tst);
    // std::cerr << "After dynamics: cice0=" << m_conc[itar] << " cice=" << ciceUpd[itar] << "
    // hice0=" << m_thick[itar] /*<< " htrue0=" << hiceTrue0[itar] */<< " htrue=" <<
    // hiceTrueUpd[itar] << std::endl;

    updatePrognosticFields();
    // std::cerr << "After updProF: cice0=" << m_conc[itar] << " cice=" << ciceUpd[itar] << "
    // hice0=" << m_thick[itar] /*<< " htrue0=" << hiceTrue0[itar] */<< " htrue=" <<
    // hiceTrueUpd[itar] << std::endl;

    // Take the updated values of the true ice and snow thicknesses, and reset hice0 and hsnow0
    // IceGrowth updates its own fields during update
    iceGrowth.update(tst);
    // std::cerr << "After IGupd: cice0=" << m_conc[itar] << " cice=" << ciceUpd[itar] << " hice0="
    // << m_thick[itar] /*<< " htrue0=" << hiceTrue0[itar] */<< " htrue=" << hiceTrueUpd[itar] <<
    // std::endl;
    updatePrognosticFields();
    // std::cerr << "After updProF2: cice0=" << m_conc[itar] << " cice=" << ciceUpd[itar] << "
    // hice0=" << m_thick[itar] /*<< " htrue0=" << hiceTrue0[itar] */<< " htrue=" <<
    // hiceTrueUpd[itar] << std::endl;

    pOcnBdy->updateAfter(tst);
}

void PrognosticData::updatePrognosticFields()
{
    ModelArrayRef<SharedArray::H_ICE, MARBackingStore, RO> hiceTrueUpd(getSharedArray());
    ModelArrayRef<SharedArray::C_ICE, MARBackingStore, RO> ciceUpd(getSharedArray());
    ModelArrayRef<SharedArray::H_SNOW, MARBackingStore, RO> hsnowTrueUpd(getSharedArray());
    ModelArrayRef<SharedArray::T_ICE, MARBackingStore, RO> ticeUpd(getSharedArray());

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
                 { "sst", mask(*getProtectedArray().at(static_cast<size_t>(ProtectedArray::SST))) },
                 { "sss", mask(*getProtectedArray().at(static_cast<size_t>(ProtectedArray::SSS))) },
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
