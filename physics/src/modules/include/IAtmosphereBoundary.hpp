/*!
 *
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#include "include/CheckingModelComponent.hpp"
#include "include/ModelArrayAccessor.hpp"
#include "include/Time.hpp"

#ifndef IATMOSPHEREBOUNDARY_HPP
#define IATMOSPHEREBOUNDARY_HPP

namespace Nextsim {

//! An interface class for the atmospheric inputs into the ice physics.
class IAtmosphereBoundary : public CheckingModelComponent {
public:
    IAtmosphereBoundary()
        : qiaAccessor(getStore(), RW, ModelArray::Type::H, std::pair(-1e4, 1e4))
        , dqia_dtAccessor(getStore(), RW, ModelArray::Type::H)
        , qowAccessor(getStore(), RW, ModelArray::Type::H, std::pair(-1e4, 1e4))
        , sublAccessor(getStore(), RW, ModelArray::Type::H)
        , snowAccessor(getStore(), RO, ModelArray::Type::H, std::pair(0.0, 1e-3))
        , rainAccessor(getStore(), RO, ModelArray::Type::H, std::pair(0.0, 1e-3))
        , evapAccessor(getStore(), RW, ModelArray::Type::H, std::pair(-1e-3, 1e-3))
        , uwindAccessor(getStore(), RO, ModelArray::Type::U, std::pair(-100.0, 100.0))
        , vwindAccessor(getStore(), RO, ModelArray::Type::V, std::pair(-100.0, 100.0))
        , penSWAccessor(getStore(), RW, ModelArray::Type::H)
        , tauXOWAccessor(getStore(), RW, ModelArray::Type::H, std::pair(-10.0, 10.0))
        , tauYOWAccessor(getStore(), RW, ModelArray::Type::H, std::pair(-10.0, 10.0))
    {
    }
    virtual ~IAtmosphereBoundary() = default;

    std::string getName() const override { return "IAtmosphereBoundary"; }
    void setData(const ModelState::DataMap& ms) override
    {
        HField& qia = qiaAccessor.getHostRW();
        qia.resize();
        HField& dqia_dt = dqia_dtAccessor.getHostRW();
        dqia_dt.resize();
        HField& qow = qowAccessor.getHostRW();
        qow.resize();
        HField& subl = sublAccessor.getHostRW();
        subl.resize();
        HField& snow = snowAccessor.getHostRW();
        snow.resize();
        HField& rain = rainAccessor.getHostRW();
        rain.resize();
        HField& evap = evapAccessor.getHostRW();
        evap.resize();
        UField& uwind = uwindAccessor.getHostRW();
        uwind.resize();
        VField& vwind = vwindAccessor.getHostRW();
        vwind.resize();
        HField& penSW = penSWAccessor.getHostRW();
        penSW.resize();
        HField& tauXOW = tauXOWAccessor.getHostRW();
        tauXOW.resize();
        HField& tauYOW = tauYOWAccessor.getHostRW();
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
    ModelArrayStore& couplingArrays() { return m_couplingArrays; }

    ModelArrayStore m_couplingArrays;

    ModelArrayAccessor<Shared::Q_IA, RW> qiaAccessor;
    ModelArrayAccessor<Shared::DQIA_DT, RW> dqia_dtAccessor;
    ModelArrayAccessor<Shared::Q_OW, RW> qowAccessor;
    ModelArrayAccessor<Shared::SUBLIM, RW> sublAccessor;
    ModelArrayAccessor<Protected::SNOW, RW> snowAccessor;
    ModelArrayAccessor<Shared::RAIN, RW> rainAccessor;
    ModelArrayAccessor<Shared::EVAP, RW> evapAccessor;
    ModelArrayAccessor<Protected::WIND_U, RW> uwindAccessor;
    ModelArrayAccessor<Protected::WIND_V, RW> vwindAccessor;
    ModelArrayAccessor<Shared::Q_PEN_SW, RW> penSWAccessor;
    ModelArrayAccessor<Shared::OW_STRESS_X, RW>
        tauXOWAccessor; // x(east)-ward open ocean stress, Pa
    ModelArrayAccessor<Shared::OW_STRESS_Y, RW>
        tauYOWAccessor; // y(north)-ward open ocean stress, Pa
};

} // namespace Nextsim

#endif /* IATMOSPHEREBOUNDARY_HPP */
