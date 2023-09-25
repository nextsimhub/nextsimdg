/*!
 * @file IAtmosphereBoundary.hpp
 *
 * @date Sep 22, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ModelArrayRef.hpp"
#include "include/ModelComponent.hpp"
#include "include/Time.hpp"

#ifndef IATMOSPHEREBOUNDARY_HPP
#define IATMOSPHEREBOUNDARY_HPP

namespace Nextsim {

namespace CouplingFields {
constexpr TextTag SUBL = "SUBL"; // sublimation mass flux kg s⁻¹ m⁻²
constexpr TextTag SNOW = "SNOW"; // snowfall mass flux kg s⁻¹ m⁻²
constexpr TextTag RAIN = "RAIN"; // rainfall mass flux kg s⁻¹ m⁻²
constexpr TextTag EVAP = "EVAP"; // evaporation mass flux kg s⁻¹ m⁻²
constexpr TextTag WIND_U = "WIND_U"; // x-aligned wind component m s⁻¹
constexpr TextTag WIND_V = "WIND_V"; // y-aligned wind component m s⁻¹

}
//! An interface class for the atmospheric inputs into the ice physics.
class IAtmosphereBoundary : public ModelComponent {
public:
    IAtmosphereBoundary()
        : qia(ModelArray::Type::H)
        , dqia_dt(ModelArray::Type::H)
    , qow(ModelArray::Type::H)
    , subl(ModelArray::Type::H)
    , snow(ModelArray::Type::H)
    , rain(ModelArray::Type::H)
    , evap(ModelArray::Type::H)
    , emp(ModelArray::Type::H)
    , uwind(ModelArray::Type::U)
    , vwind(ModelArray::Type::V)
    {
        m_couplingArrays.registerArray(CouplingFields::SUBL, &subl, RW);
        m_couplingArrays.registerArray(CouplingFields::SNOW, &snow, RW);
        m_couplingArrays.registerArray(CouplingFields::RAIN, &rain, RW);
        m_couplingArrays.registerArray(CouplingFields::EVAP, &evap, RW);
        m_couplingArrays.registerArray(CouplingFields::WIND_U, &uwind, RW);
        m_couplingArrays.registerArray(CouplingFields::WIND_V, &vwind, RW);

        getStore().registerArray(Shared::Q_IA, &qia, RW);
        getStore().registerArray(Shared::DQIA_DT, &dqia_dt, RW);
        getStore().registerArray(Shared::Q_OW, &qow, RW);
        getStore().registerArray(Shared::SUBLIM, &subl, RW);
        getStore().registerArray(Protected::SNOW, &snow, RO);
        getStore().registerArray(Shared::Q_PEN_SW, &penSW, RW);
        getStore().registerArray(Protected::EVAP_MINUS_PRECIP, &emp, RO);
        getStore().registerArray(Protected::WIND_U, &uwind, RO);
        getStore().registerArray(Protected::WIND_V, &vwind, RO);
    }
    virtual ~IAtmosphereBoundary() = default;

    ModelState getState() const override { return ModelState(); }
    ModelState getState(const OutputLevel&) const override { return getState(); }

    std::string getName() const override { return "IAtmosphereBoundary"; }
    void setData(const ModelState::DataMap& ms) override
    {
        qia.resize();
        dqia_dt.resize();
        qow.resize();
        subl.resize();
        snow.resize();
        rain.resize();
        evap.resize();
        emp.resize();
        uwind.resize();
        vwind.resize();
        penSW.resize();
    }
    virtual void update(const TimestepTime& tst) { }

protected:

    ModelArrayReferenceStore& couplingArrays() { return m_couplingArrays; }

    HField qia;
    HField dqia_dt;
    HField qow;
    HField subl;
    HField snow;
    HField rain;
    HField evap;
    HField emp;
    UField uwind;
    VField vwind;
    HField penSW;

    ModelArrayReferenceStore m_couplingArrays;
};

} // namespace Nextsim

#endif /* IATMOSPHEREBOUNDARY_HPP */
