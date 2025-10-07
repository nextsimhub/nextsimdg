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
        // Make sure the outputs have the right size
        uOut.resize();
        vOut.resize();

        uOut = uIn * cosAngle - vIn * sinAngle;
        vOut = uIn * sinAngle + vIn * cosAngle;
    };

    void rotateVectorToGreenland(
        const ModelArray& uIn, const ModelArray& vIn, ModelArray& uOut, ModelArray& vOut) const
    {
        // Make sure the outputs have the right size
        uOut.resize();
        vOut.resize();

        uOut = uIn * cosAngle + vIn * sinAngle;
        vOut = -uIn * sinAngle + vIn * cosAngle;
    };

    ModelArray cosAngle, sinAngle;

    /*!
     * Interpolate fields on U and V points of a C-grid to the centre point of the grid. The U and V
     * points should be placed on the right and above the centre point.
     *
     * @param uIn Input U field on the U-point
     * @param vIn Input V field on the V-point
     * @param uOut U field on centre point
     * @param vOut V field on centre point
     */
    void interpCURtoA(const UField& uIn, const VField& vIn, UField& uOut, VField& vOut) const
    {
        // Make sure the outputs have the right size
        uOut.resize();
        vOut.resize();

        // First the interior
        for (size_t i = 1; i < uIn.dimensions()[0]; ++i) {
            for (size_t j = 1; j < uIn.dimensions()[1]; ++j) {
                uOut(i, j) = 0.5 * (uIn(i, j) + uIn(i - 1, j));
                vOut(i, j) = 0.5 * (vIn(i, j) + vIn(i, j - 1));
            }
        }

        // Then the bottom row
        for (size_t i = 1; i < uIn.dimensions()[0]; ++i)
            uOut(i, 0) = 0.5 * (uIn(i, 0) + uIn(i - 1, 0));

        // Then the first column
        for (size_t j = 1; j < uIn.dimensions()[1]; ++j)
            vOut(0, j) = 0.5 * (vIn(0, j) + vIn(0, j - 1));

        // Then the bottom left corner
        uOut(0, 0) = 0.5 * uIn(0, 0);
        vOut(0, 0) = 0.5 * vIn(0, 0);
    }

    /*!
     * Interpolate fields on the centre point of the grid to the U and V points of a C-grid. The U
     * and V points should be placed on the right and above the centre point.
     *
     * @param uIn Input U field on the centre point
     * @param vIn Input V field on the centre point
     * @param uOut U field on the U-point
     * @param vOut V field on the V-point
     */
    void interpAtoCUR(const UField& uIn, const VField& vIn, UField& uOut, VField& vOut) const
    {
        // Make sure the outputs have the right size
        uOut.resize();
        vOut.resize();

        // First the interior
        for (size_t i = 1; i < uIn.dimensions()[0] - 1; ++i) {
            for (size_t j = 1; j < uIn.dimensions()[1] - 1; ++j) {
                uOut(i, j) = 0.5 * (uIn(i, j) + uIn(i + 1, j));
                vOut(i, j) = 0.5 * (vIn(i, j) + vIn(i, j + 1));
            }
        }

        // Then the top row
        for (size_t i = 1; i < uIn.dimensions()[0]; ++i)
            uOut(i, 0) = 0.5 * (uIn(i, 0) + uIn(i - 1, 0));

        // Then the first column
        for (size_t j = 1; j < uIn.dimensions()[1]; ++j)
            vOut(0, j) = 0.5 * (vIn(0, j) + vIn(0, j - 1));

        // Then the bottom left corner
        uOut(0, 0) = 0.5 * uIn(0, 0);
        vOut(0, 0) = 0.5 * vIn(0, 0);
    }
};

}

#endif // OASISCOUPLED_HPP
