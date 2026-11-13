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
        , qowAccessor(getStore(), RW, ModelArray::Type::H)
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
        qia.reinitialize();
        HField& dqia_dt = dqia_dtAccessor.getHostRW();
        dqia_dt.reinitialize();
        HField& qow = qowAccessor.getHostRW();
        qow.reinitialize();
        HField& subl = sublAccessor.getHostRW();
        subl.reinitialize();
        HField& snow = snowAccessor.getHostRW();
        snow.reinitialize();
        HField& rain = rainAccessor.getHostRW();
        rain.reinitialize();
        HField& evap = evapAccessor.getHostRW();
        evap.reinitialize();
        UField& uwind = uwindAccessor.getHostRW();
        uwind.reinitialize();
        VField& vwind = vwindAccessor.getHostRW();
        vwind.reinitialize();
        HField& penSW = penSWAccessor.getHostRW();
        penSW.reinitialize();
        HField& tauXOW = tauXOWAccessor.getHostRW();
        tauXOW.reinitialize();
        HField& tauYOW = tauYOWAccessor.getHostRW();
        tauYOW.reinitialize();

        addChecks({
            { "qia", qiaAccessor },
            { "qow", qowAccessor },
            { "snow", snowAccessor },
            { "rain", rainAccessor },
            { "evap", evapAccessor },
            { "uwind", uwindAccessor },
            { "vwind", vwindAccessor },
            { "tauXOW", tauXOWAccessor },
            { "tauYOW", tauYOWAccessor },
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
