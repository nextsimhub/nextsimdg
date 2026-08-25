/*!
 *
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IFLUXCALCULATION_HPP
#define IFLUXCALCULATION_HPP

#include "include/Configured.hpp"
#include "include/ModelArrayAccessor.hpp"
#include "include/ModelComponent.hpp"
#include "include/ModelState.hpp"

namespace Nextsim {
//! An interface class for calculating ocean/ice-atmosphere fluxes
class IFluxCalculation : public ModelComponent {
public:
    IFluxCalculation()
        : qowAccessor(getStore())
        , evapAccessor(getStore())
        , sublAccessor(getStore())
        , qiaAccessor(getStore())
        , penSWAccessor(getStore())
        , dqia_dtAccessor(getStore())
        , Q_sw_owAccessor(getStore())
        , Q_sw_baseAccessor(getStore())
        , tau_x_owAccessor(getStore())
        , tau_y_owAccessor(getStore())
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
    ModelArrayAccessor<Shared::Q_OW, RW> qowAccessor;
    ModelArrayAccessor<Shared::EVAP, RW> evapAccessor;
    ModelArrayAccessor<Shared::SUBLIM, RW> sublAccessor;
    ModelArrayAccessor<Shared::Q_IA, RW> qiaAccessor;
    ModelArrayAccessor<Shared::Q_PEN_SW, RW> penSWAccessor;
    ModelArrayAccessor<Shared::DQIA_DT, RW> dqia_dtAccessor;
    ModelArrayAccessor<Shared::Q_SW_OW, RW> Q_sw_owAccessor;
    ModelArrayAccessor<Shared::Q_SW_BASE, RW> Q_sw_baseAccessor;
    ModelArrayAccessor<Shared::OW_STRESS_X, RW> tau_x_owAccessor;
    ModelArrayAccessor<Shared::OW_STRESS_Y, RW> tau_y_owAccessor;
};
}
#endif /* IFLUXCALCULATION_HPP */
