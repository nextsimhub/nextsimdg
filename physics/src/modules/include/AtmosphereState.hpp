/*!
 * @file AtmosphereState.hpp
 *
 * @date May 9, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef ATMOSPHERESTATE_HPP
#define ATMOSPHERESTATE_HPP

#include "include/ModelComponent.hpp"
#include "include/Time.hpp"

namespace Nextsim {

//! A class providing an interface to the atmosphere climatologies or coupled model.
class AtmosphereState : public ModelComponent {
public:
    AtmosphereState();
    virtual ~AtmosphereState() = default;

    void setData(const ModelState::DataMap&) override;
    ModelState getState() const override;
    ModelState getState(const OutputLevel&) const override;
    ModelState getStateRecursive(const OutputSpec& os) const override;

    std::string getName() const override;
    std::unordered_set<std::string> hFields() const override;

    /*!
     * @brief Updates the atmosphere state.
     *
     * @details Performs any common calculations, then any implementation
     * specific updates.
     *
     * @param tStep The object containing the timestep start and duration times.
     */
    void update(const TimestepTime& tStep);

protected:
    HField tair;
    HField tdew;
    HField pair;
    HField rmix;
    HField sw_in;
    HField lw_in;
    HField snowfall;
    HField evap_minus_precip;

    HField windSpeed;

    virtual void updateSpecial(const TimestepTime&) = 0;
};

} /* namespace Nextsim */

#endif /* ATMOSPHERESTATE_HPP */
