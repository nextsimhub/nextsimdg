/*!
 * @file IAtmosphereBoundary.hpp
 *
 * @date Sep 12, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IATMOSPHEREBOUNDARY_HPP
#define IATMOSPHEREBOUNDARY_HPP

#include "include/ModelComponent.hpp"

namespace Nextsim {

class IAtmosphereBoundary : public ModelComponent {
    IAtmosphereBoundary()
    {
        registerSharedArray(SharedArray::Q_IA, &qia);
        registerSharedArray(SharedArray::DQIA_DT, &dqia_dt);
        registerProtectedArray(ProtectedArray::SW_IN, &qsw);
        registerProtectedArray(ProtectedArray::LW_IN, &qlw);
        registerSharedArray(SharedArray::SUBLIM, &subl);

    }
    virtual ~IAtmosphereBoundary() = default;

    void setData(const ModelState::DataMap&) override;
    ModelState getState() const override;
    ModelState getState(const OutputLevel&) const override;
    ModelState getStateRecursive(const OutputSpec& os) const override;

    std::string getName() const override { return "IAtmosphereBoundary"; }
    std::unordered_set<std::string> hFields() const override;

protected:
    HField qia; // Ice-atmosphere heat flux, W m⁻²
    HField dqia_dt; // Derivative of qia w.r.t. ice surface temperature
    HField qsw; // Solar shortwave radiation, W m⁻²
    HField qlw; // Non-Solar radiative heat flux, W m⁻²
    HField subl; // Ice sublimative mass flux, kg m⁻²
    HField precip; // total precipitation flux, kg m⁻²
    UField u; // x(east)-ward wind, m s⁻¹
    VField v; // y(north)-ward wind, m s⁻¹
};

} /* namespace Nextsim */

#endif /* IATMOSPHEREBOUNDARY_HPP */
