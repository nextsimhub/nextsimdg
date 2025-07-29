/*!
 *
 * @author  Tim Spain <timothy.spain@nersc.no>
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
        , qow(ModelArray::Type::H, { -1e4, 1e4 })
        , subl(ModelArray::Type::H)
        , snow(ModelArray::Type::H, { 0, 1e-3 })
        , rain(ModelArray::Type::H, { 0, 1e-3 })
        , evap(ModelArray::Type::H, { -1e-3, 1e-3 })
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
        getStore().registerArray(Shared::EVAP, &evap, RW);
        getStore().registerArray(Shared::RAIN, &rain, RO);
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
        uwind.resize();
        vwind.resize();
        penSW.resize();
        tauXOW.resize();
        tauYOW.resize();

        addChecks({
            { "qia", &qia },
            { "qow", &qow },
            { "snow", &snow },
            { "rain", &rain },
            { "evap", &evap },
            { "uwind", &uwind },
            { "vwind", &vwind },
            { "tauXOW", &tauXOW },
            { "tauYOW", &tauYOW },
        });
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
    UField uwind;
    VField vwind;
    HField penSW;
    HField tauXOW; // x(east)-ward open ocean stress, Pa
    HField tauYOW; // y(north)-ward open ocean stress, Pa

    ModelArrayReferenceStore m_couplingArrays;
};

} // namespace Nextsim

#endif /* IATMOSPHEREBOUNDARY_HPP */
