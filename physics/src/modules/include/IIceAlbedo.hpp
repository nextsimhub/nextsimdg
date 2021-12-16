/*
 * @file IIceAlbedo.hpp
 *
 * @date Sep 28, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_IICEALBEDO_HPP
#define SRC_INCLUDE_IICEALBEDO_HPP

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
    virtual double albedo(double temperature, double snowThickness) = 0;
};
}
#endif /* SRC_INCLUDE_IICEALBEDO_HPP */
