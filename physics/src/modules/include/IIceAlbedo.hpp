/*
 * @file IIceAlbedo.hpp
 *
 * @date Sep 28, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IICEALBEDO_HPP
#define IICEALBEDO_HPP

#include "include/Time.hpp"

namespace Nextsim {
//! The interface class for ice albedo calculation.
class IIceAlbedo {
public:
    virtual ~IIceAlbedo() = default;
    /*!
     * @brief A virtual function that calculates the ice surface short wave
     * albedo.
     *
     * @param temperature The temperature of the ice surface.
     * @param snowThickness The true snow thickness on top of the ice.
     */
    virtual std::tuple<double, double> albedo(double temperature, double snowThickness, double i0)
        = 0;

    /*!
     * Sets the time parameter for the implementation, if it is time dependent.
     * @param time The desired TimePoint.
     */
    virtual void setTime(const TimePoint& tp) { }
};
}
#endif /* IICEALBEDO_HPP */
