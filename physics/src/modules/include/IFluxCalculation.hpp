/*!
 * @file IFluxCalculation.hpp
 *
 * @date Apr 29, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IFLUXCALCULATION_HPP
#define IFLUXCALCULATION_HPP

#include "include/ModelComponent.hpp"
#include "include/ModelState.hpp"

namespace Nextsim {
class IFluxCalculation : public ModelComponent {
public:
    IFluxCalculation()
    {
        // register shared arrays
        registerSharedArray(SharedArray::Q_OW, &qow);
        registerSharedArray(SharedArray::SUBLIM, &subl);
        registerSharedArray(SharedArray::Q_IA, &qia);
        registerSharedArray(SharedArray::DQIA_DT, &dqia_dt);
        registerSharedArray(SharedArray::Q_IO, &qio);
    }
    virtual ~IFluxCalculation() = default;

    void setData(const ModelState&) override { }

    ModelState getState() const override { return ModelState(); }
    ModelState getState(const OutputLevel&) const override { return getState(); }

    std::string getName() const override { return "IFluxCalculation"; }

    std::set<std::string> hFields() const override
    {
        return { "qow", "subl", "qia", "dqia_dt", "qio" };
    }

    virtual void update(const TimestepTime&) = 0;

protected:
    // All fluxes are positive upwards, including incident radiation fluxes
    // Owned, shared fields
    HField qow; // Open water heat flux [W m⁻²]
    HField subl; // Ice sublimative mass flux [kg m⁻²]
    HField qia; // Ice-atmosphere heat flux [W m⁻²]
    HField dqia_dt; // Derivative of qia w.r.t. ice surface temperature
    HField qio; // Ice-ocean heat flux [W m⁻²]
};
}
#endif /* IFLUXCALCULATION_HPP */
