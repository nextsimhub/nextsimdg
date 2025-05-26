/*!
 * @file IFluxCalculation.hpp
 *
 * @date 23 May 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IFLUXCALCULATION_HPP
#define IFLUXCALCULATION_HPP

#include "include/Configured.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/ModelComponent.hpp"
#include "include/ModelState.hpp"

namespace Nextsim {
//! An interface class for calculating ocean/ice-atmosphere fluxes
class IFluxCalculation : public ModelComponent {
public:
    IFluxCalculation()
        : qow(getStore())
        , evap(getStore())
        , subl(getStore())
        , qia(getStore())
        , penSW(getStore())
        , dqia_dt(getStore())
        , Q_sw_ow(getStore())
        , Q_sw_base(getStore())
        , tau_x_ow(getStore())
        , tau_y_ow(getStore())
    {
    }
    virtual ~IFluxCalculation() = default;

    void setData(const ModelState::DataMap& ms) override { }

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
    ModelArrayRef<Shared::Q_OW, RW> qow;
    ModelArrayRef<Shared::EVAP, RW> evap;
    ModelArrayRef<Shared::SUBLIM, RW> subl;
    ModelArrayRef<Shared::Q_IA, RW> qia;
    ModelArrayRef<Shared::Q_PEN_SW, RW> penSW;
    ModelArrayRef<Shared::DQIA_DT, RW> dqia_dt;
    ModelArrayRef<Shared::Q_SW_OW, RW> Q_sw_ow;
    ModelArrayRef<Shared::Q_SW_BASE, RW> Q_sw_base;
    ModelArrayRef<Shared::OW_STRESS_X, RW> tau_x_ow;
    ModelArrayRef<Shared::OW_STRESS_Y, RW> tau_y_ow;
};
}
#endif /* IFLUXCALCULATION_HPP */
