/*!
 * @file ILateralIceSpread.hpp
 *
 * @date Apr 5, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef ILATERALICESPREAD_HPP
#define ILATERALICESPREAD_HPP

#include "include/ModelComponent.hpp"
#include "include/Time.hpp"

namespace Nextsim {

class ILateralIceSpread : public ModelComponent {
public:
    virtual ~ILateralIceSpread() = default;

    std::string getName() const override { return "LateralIceSpread"; }
    void setData(const ModelState& ms) override { }
    ModelState getState() const override { return ModelState(); }
    ModelState getState(const OutputLevel&) const override { return getState(); }
    virtual void freeze(const TimestepTime& tstep, double hice, double hsnow, double deltaHi,
        double newIce, double& cice, double& qow, double& deltaCfreeze)
        = 0;
    virtual void melt(const TimestepTime& tstep, double hice, double hsnow, double deltaHi,
        double& cice, double& qow, double& deltaCmelt)
        = 0;

protected:
    ILateralIceSpread()
    {
        registerModule();

        ModelComponent::requestSharedArray(SharedArray::C_ICE, &cice);
        ModelComponent::requestSharedArray(SharedArray::Q_OW, &qow);

        ModelComponent::requestProtectedArray(SharedArray::H_ICE, &hice);
        ModelComponent::requestProtectedArray(SharedArray::H_SNOW, &hsnow);
        ModelComponent::requestProtectedArray(SharedArray::DELTA_HICE, &deltaHi);

        ModelComponent::registerSharedArray(SharedArray::DELTA_CICE, &deltaCi);
    }

    HField* cice; // From IceGrowth
    HField* qow; // From FluxCalculation
    const HField* hice; // From IceGrowth
    const HField* hsnow; // From Ice Growth?
    const HField* deltaHi; // From Vertical Ice Growth

    // Owned, shared arrays
    HField deltaCi;
};

} /* namespace Nextsim */

#endif /* ILATERALICESPREAD_HPP */
