/*!
 * @file IHorizontalIceSpread.hpp
 *
 * @date Apr 5, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IHORIZONTALICESPREAD_HPP
#define IHORIZONTALICESPREAD_HPP

#include "include/ModelComponent.hpp"
#include "include/Time.hpp"

namespace Nextsim {

class IHorizontalIceSpread : public ModelComponent {
public:
    virtual ~IHorizontalIceSpread();

    std::string getName() const override {return "HorizontalIceSpread"; }
    void setData(const ModelState& ms) override { }
    ModelState getState() const override { return ModelState(); }
    ModelState getState(const OutputLevel&) const override { return getState(); }
    virtual void update(const TimestepTime&) = 0;
protected:
    IHorizontalIceSpread()
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

#endif /* IHORIZONTALICESPREAD_HPP */
