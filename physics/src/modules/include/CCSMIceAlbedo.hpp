/*!
 * @file CCSMIceAlbedo.hpp
 *
 * @date Sep 28, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_CCSMICEALBEDO_HPP
#define SRC_INCLUDE_CCSMICEALBEDO_HPP

#include "include/Configured.hpp"
#include "IIceAlbedo.hpp"

namespace Nextsim {

//! The implementation class for the CCSM calculation of ice surface albedo.
class CCSMIceAlbedo : public IIceAlbedo, public Configured<CCSMIceAlbedo> {
public:
    /*!
     * @brief Calculates the CCSM ice surface short wave albedo.
     *
     * @param temperature The temperature of the ice surface.
     * @param snowThickness The true snow thickness on top of the ice.
     */
    double albedo(double temperature, double snowThickness) override;

    void configure() override;

private:
    static double iceAlbedo;
    static double snowAlbedo;
};

}

#endif /* SRC_INCLUDE_CCSMICEALBEDO_HPP */
