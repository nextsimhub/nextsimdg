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
class IAtmosphereBoundary : public ModelComponent {
public:
    IAtmosphereBoundary()
    {
        m_couplingArrays.resize(static_cast<size_t>(CouplingFields::COUNT));
        m_couplingArrays[static_cast<size_t>(CouplingFields::SUBL)] = &subl;
        m_couplingArrays[static_cast<size_t>(CouplingFields::PRECIP)] = &precip;
        m_couplingArrays[static_cast<size_t>(CouplingFields::EVAP)] = &evap;
        m_couplingArrays[static_cast<size_t>(CouplingFields::SW_IN)] = &sw_in;
        m_couplingArrays[static_cast<size_t>(CouplingFields::LW_IN)] = &lw_in;
        m_couplingArrays[static_cast<size_t>(CouplingFields::WIND_U)] = &uwind;
        m_couplingArrays[static_cast<size_t>(CouplingFields::WIND_V)] = &vwind;

        registerSharedArray(SharedArray::Q_IA, &qia);
        registerSharedArray(SharedArray::DQIA_DT, &dqia_dt);
        registerSharedArray(SharedArray::SUBLIM, &subl);
        registerProtectedArray(ProtectedArray::SW_IN, &sw_in);
        registerProtectedArray(ProtectedArray::LW_IN, &lw_in);
    }
    virtual ~IAtmosphereBoundary() = default;

    ModelState getState() const override { return ModelState(); }
    ModelState getState(const OutputLevel&) const override { return getState(); }

    std::string getName() const override { return "IAtmosphereBoundary"; }
    virtual void setData(const ModelState::DataMap& ms)
    {
        qia.resize();
        dqia_dt.resize();
        subl.resize();
        precip.resize();
        evap.resize();
        sw_in.resize();
        lw_in.resize();
        uwind.resize();
        vwind.resize();
    }
    virtual void update(const TimestepTime& tst) { }
    protected:
        enum class CouplingFields {
            SUBL, // sublimation mass flux kg s⁻¹ m⁻²
            PRECIP, // precipitation mass flux kg s⁻¹ m⁻²
            EVAP, // evaporation mass flux kg s⁻¹ m⁻²
            SW_IN, // solar radiation W m⁻²
            LW_IN, // All non-solar radiation W m⁻²
            WIND_U, // x-aligned wind component m s⁻¹
            WIND_V, // y-aligned wind component m s⁻¹
            COUNT
        };

        const MARBackingStore& couplingArrays() { return m_couplingArrays; }

        HField qia;
        HField dqia_dt;
        HField subl;
        HField precip;
        HField sw_in;
        HField lw_in;
        HField evap;
        UField uwind;
        VField vwind;

        MARBackingStore m_couplingArrays;
};

} // namespace Nextsim

#endif /* IATMOSPHEREBOUNDARY_HPP */
