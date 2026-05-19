/*!
 *
 * @author  Einar Ólason <einar.olason@nersc.no>
 */

#ifndef IDAMAGEHEALING_HPP
#define IDAMAGEHEALING_HPP

#include "include/ModelArrayAccessor.hpp"
#include "include/ModelComponent.hpp"

namespace Nextsim {
//! An interface class for modules controlling damage healing
class IDamageHealing : public ModelComponent {
public:
    virtual ~IDamageHealing() = default;

    std::string getName() const override { return "DamageHealing"; }
    void setData(const ModelState::DataMap& ms) override { }
    /*!
     * Updates the ice damage based on lateral ice growth and healing
     *
     * @param tStep The object containing the timestep start and duration times.
     */
    virtual void update(const TimestepTime& tsTime) = 0;

protected:
    IDamageHealing()
        : ciceAccessor(getStore())
        , deltaCiAccessor(getStore())
        , damageAccessor(getStore())
    {
    }

    ModelArrayAccessor<Shared::C_ICE_DG, RO> ciceAccessor; // From prognostic data
    ModelArrayAccessor<Shared::DELTA_CICE, RO> deltaCiAccessor; // From LateralIceSpread
    ModelArrayAccessor<Shared::DAMAGE, RW> damageAccessor; // From prognostic data
};

} /* namespace Nextsim */

#endif // IDAMAGEHEALING_HPP
