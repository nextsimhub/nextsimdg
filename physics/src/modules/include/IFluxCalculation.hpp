/*!
 * @file IFluxCalculation.hpp
 *
 * @date Apr 29, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IFLUXCALCULATION_HPP
#define IFLUXCALCULATION_HPP

#include "include/ModelArrayRef.hpp"
#include "include/Configured.hpp"
#include "include/ModelComponent.hpp"
#include "include/ModelState.hpp"

namespace Nextsim {
//! An interface class for calculating ocean/ice-atmosphere fluxes
class IFluxCalculation : public ModelComponent {
public:
    IFluxCalculation()
        : qow(getStore())
        , subl(getStore())
        , qia(getStore())
        , penSW(getStore())
        , dqia_dt(getStore())
    {
    }
    virtual ~IFluxCalculation() = default;

    void setData(const ModelState::DataMap& ms) override { }

    ModelState getState() const override { return ModelState(); }
    ModelState getState(const OutputLevel&) const override { return getState(); }
    ModelState getStateRecursive(const OutputSpec& os) const override
    {
        return os ? getState() : ModelState();
    }

    std::string getName() const override { return "IFluxCalculation"; }

    /*!
     * Updates the flux calculation for the timestep.
     *
     * @param tStep The object containing the timestep start and duration times.
     */
    virtual void update(const TimestepTime&) = 0;

protected:
    // All fluxes are positive upwards, including incident radiation fluxes
    // The flux fields are owned by IAtmosphereBoundary
    ModelArrayRef<Shared::Q_PEN_SW, RW> penSW;
    ModelArrayRef<Shared::Q_OW, RW> qow;
    ModelArrayRef<Shared::SUBLIM, RW> subl;
    ModelArrayRef<Shared::Q_IA, RW> qia;
    ModelArrayRef<Shared::DQIA_DT, RW> dqia_dt;
};
}
#endif /* IFLUXCALCULATION_HPP */
