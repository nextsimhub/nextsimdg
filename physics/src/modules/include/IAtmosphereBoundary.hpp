/*!
 * @file IAtmosphereBoundary.hpp
 *
 * @date Sep 22, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "../../../../core/src/include/ModelArrayRef.hpp"
#include "include/ModelComponent.hpp"
#include "include/Time.hpp"

#ifndef IATMOSPHEREBOUNDARY_HPP
#define IATMOSPHEREBOUNDARY_HPP

namespace Nextsim {

//! An interface class for the atmospheric inputs into the ice physics.
class IAtmosphereBoundary : public ModelComponent {
public:
    IAtmosphereBoundary()
    {
        m_couplingArrays.resize(static_cast<size_t>(CouplingFields::COUNT));
        m_couplingArrays[static_cast<size_t>(CouplingFields::SUBL)] = &subl;
        m_couplingArrays[static_cast<size_t>(CouplingFields::SNOW)] = &snow;
        m_couplingArrays[static_cast<size_t>(CouplingFields::RAIN)] = &rain;
        m_couplingArrays[static_cast<size_t>(CouplingFields::EVAP)] = &evap;
        m_couplingArrays[static_cast<size_t>(CouplingFields::WIND_U)] = &uwind;
        m_couplingArrays[static_cast<size_t>(CouplingFields::WIND_V)] = &vwind;

        getStore().registerArray(Shared::Q_IA, &qia, RW);
        getStore().registerArray(Shared::DQIA_DT, &dqia_dt, RW);
        getStore().registerArray(Shared::Q_OW, &qow, RW);
        getStore().registerArray(Shared::SUBLIM, &subl, RW);
        getStore().registerArray(Protected::SNOW, &snow);
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
        uwind.resize();
        vwind.resize();
    }
    virtual void update(const TimestepTime& tst) { }

protected:
    enum class CouplingFields {
        SUBL, // sublimation mass flux kg s⁻¹ m⁻²
        SNOW, // snowfall mass flux kg s⁻¹ m⁻²
        RAIN, // rainfall mass flux kg s⁻¹ m⁻²
        EVAP, // evaporation mass flux kg s⁻¹ m⁻²
        WIND_U, // x-aligned wind component m s⁻¹
        WIND_V, // y-aligned wind component m s⁻¹
        COUNT
    };

    const MARBackingStore& couplingArrays() { return m_couplingArrays; }

    HField qia;
    HField dqia_dt;
    HField qow;
    HField subl;
    HField snow;
    HField rain;
    HField evap;
    UField uwind;
    VField vwind;

    MARBackingStore m_couplingArrays;
};

} // namespace Nextsim

#endif /* IATMOSPHEREBOUNDARY_HPP */
