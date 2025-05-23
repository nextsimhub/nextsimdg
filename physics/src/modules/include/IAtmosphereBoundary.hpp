/*!
 * @file IAtmosphereBoundary.hpp
 *
 * @date 23 May 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/CheckingModelComponent.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/Time.hpp"

#ifndef IATMOSPHEREBOUNDARY_HPP
#define IATMOSPHEREBOUNDARY_HPP

namespace Nextsim {

//! An interface class for the atmospheric inputs into the ice physics.
class IAtmosphereBoundary : public CheckingModelComponent {
public:
    IAtmosphereBoundary()
        : qia(ModelArray::Type::H, { -1e4, 1e4 })
        , dqia_dt(ModelArray::Type::H)
        , qow(ModelArray::Type::H, { -1e3, 1e3 })
        , subl(ModelArray::Type::H)
        , snow(ModelArray::Type::H, { 0, 1e-3 })
        , rain(ModelArray::Type::H, { 0, 1e-3 })
        , evap(ModelArray::Type::H, { 0, 1e-3 })
        , emp(ModelArray::Type::H, { -1e-3, 1e-3 })
        , uwind(ModelArray::Type::U, { -100, 100 })
        , vwind(ModelArray::Type::V, { -100, 100 })
        , penSW(ModelArray::Type::H)
        , tauXOW(ModelArray::Type::H, { -10, 10 })
        , tauYOW(ModelArray::Type::H, { -10, 10 })
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
        getStore().registerArray(Shared::OW_STRESS_X, &tauXOW, RW);
        getStore().registerArray(Shared::OW_STRESS_Y, &tauYOW, RW);
        getStore().registerArray(Protected::SNOW, &snow, RO);
        getStore().registerArray(Protected::EVAP_MINUS_PRECIP, &emp, RO);
        getStore().registerArray(Protected::WIND_U, &uwind, RO);
        getStore().registerArray(Protected::WIND_V, &vwind, RO);
        getStore().registerArray(Shared::Q_PEN_SW, &penSW, RW);
    }
    virtual ~IAtmosphereBoundary() = default;

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
        tauXOW.resize();
        tauYOW.resize();

        fieldsToCheck.emplace_back("qia", &qia);
        fieldsToCheck.emplace_back("qow", &qow);
        fieldsToCheck.emplace_back("snow", &snow);
        fieldsToCheck.emplace_back("rain", &rain);
        fieldsToCheck.emplace_back("evap", &evap);
        fieldsToCheck.emplace_back("emp", &emp);
        fieldsToCheck.emplace_back("uwind", &uwind);
        fieldsToCheck.emplace_back("vwind", &vwind);
        fieldsToCheck.emplace_back("tauXOW", &tauXOW);
        fieldsToCheck.emplace_back("tauYOW", &tauYOW);
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
    HField tauXOW; // x(east)-ward open ocean stress, Pa
    HField tauYOW; // y(north)-ward open ocean stress, Pa

    ModelArrayReferenceStore m_couplingArrays;
};

} // namespace Nextsim

#endif /* IATMOSPHEREBOUNDARY_HPP */
