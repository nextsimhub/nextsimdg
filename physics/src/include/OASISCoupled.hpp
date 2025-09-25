/*!
 * @file OASISCoupled.hpp
 *
 * @date 09 Sep 2024
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#ifndef OASISCOUPLED_HPP
#define OASISCOUPLED_HPP

#ifdef USE_OASIS
#include <oasis_c.h>
#endif

#include "include/ModelArray.hpp"
#include "include/Time.hpp"

namespace Nextsim {

class OASISCoupled {
public:
    virtual std::string getName() const { return "OASISCoupled"; }

#ifdef USE_OASIS
    int OASISTime;

    // Set the "OASIS time" (seconds since start) to zero
    OASISCoupled() { OASISTime = 0; }

    // Increment the "OASIS" time by the number of seconds in the time step
    // Could be any time unit
    // Must be called at the end of the child class' update or updateAfter call.
    void updateOASISTime(const TimestepTime& tst) { OASISTime += tst.step.seconds(); }
#else
    const std::string NoOASISError = std::string(": OASIS support not compiled in.\n");
#endif

protected:
    void rotateVectorFromGreenland(
        const ModelArray& uIn, const ModelArray& vIn, ModelArray& uOut, ModelArray& vOut) const
    {
        uOut = uIn * cosAngle - vIn * sinAngle;
        vOut = uIn * sinAngle + vIn * cosAngle;
    };

    void rotateVectorToGreenland(
        const ModelArray& uIn, const ModelArray& vIn, ModelArray& uOut, ModelArray& vOut) const
    {
        uOut = uIn * cosAngle + vIn * sinAngle;
        vOut = -uIn * sinAngle + vIn * cosAngle;
    };

    ModelArray cosAngle, sinAngle;
};

}

#endif // OASISCOUPLED_HPP
