/*!
 * @file HoHealing.hpp
 *
 * This class has no corresponding implementation, just this header file
 *
 * @date 21 Nov 2024
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#ifndef NOHEALING_HPP
#define NOHEALING_HPP

#include "include/IDamageHealing.hpp"

namespace Nextsim {

//! A class implementing no healing of damage
class NoHealing : public IDamageHealing {
public:
    NoHealing()
        : IDamageHealing()
    {
    }
    virtual ~NoHealing() = default;

    ModelState getStateRecursive(const OutputSpec& os) const override { return ModelState(); };

    // Do nothing when update is called.
    void update(const TimestepTime& tstep) override { damage = oldDamage; }
};

}

#endif // NOHEALING_HPP
