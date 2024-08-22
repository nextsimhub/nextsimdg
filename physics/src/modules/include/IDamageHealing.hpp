/*!
 * @file IDamageHealing.hpp
 *
 * @date Jun 3, 2024
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#ifndef IDAMAGEHEALING_HPP
#define IDAMAGEHEALING_HPP

#include "include/ModelArrayRef.hpp"
#include "include/ModelComponent.hpp"

namespace Nextsim {
//! An interface class for modules controlling damage healing
class IDamageHealing : public ModelComponent {
public:
    virtual ~IDamageHealing() = default;

    std::string getName() const override { return "DamageHealing"; }
    void setData(const ModelState::DataMap& ms) override {
        damage.data().resize();
    }
    ModelState getState() const override { return ModelState(); }
    ModelState getState(const OutputLevel&) const override { return getState(); }
    ModelState getStateRecursive(const OutputSpec& os) const override
    {
        return os ? getState() : ModelState();
    }
    /*!
     * Updates the ice damage based on lateral ice growth and healing
     *
     * @param tStep The object containing the timestep start and duration times.
     */
    virtual void update(const TimestepTime& tsTime) = 0;

protected:
    IDamageHealing()
        : cice(getStore())
        , deltaCi(getStore())
        , damage(getStore())
    {
    }

    ModelArrayRef<Shared::C_ICE, RO> cice; // From prognostic data
    ModelArrayRef<Shared::DELTA_CICE, RO> deltaCi; // From LateralIceSpread
    ModelArrayRef<Shared::DAMAGE, RW> damage; // From prognostic data
};

} /* namespace Nextsim */

#endif // IDAMAGEHEALING_HPP
