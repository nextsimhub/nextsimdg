/*!
 * @file SMU2IceAlbedo.hpp
 *
 * @date Sep 28, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_SMU2ICEALBEDO_HPP
#define SRC_INCLUDE_SMU2ICEALBEDO_HPP

#include "include/IIceAlbedo.hpp"

namespace Nextsim {

//! The implementation class for the SMU calculation of ice surface albedo
// with variable snow albedo.
class SMU2IceAlbedo : public IIceAlbedo {
    /*!
     * @brief Calculates the SMU ice surface short wave albedo with constant
     * snow albedo.
     *
     * @param temperature The temperature of the ice surface.
     * @param snowThickness The true snow thickness on top of the ice.
     * @param i0 The fraction of short-wave radiation that can penetrate bare ice (not taking snow
     * cover into account).
     */
    std::tuple<double, double> surfaceShortWaveBalance(double temperature, double snowThickness, double i0);
};

}

#endif /* SRC_INCLUDE_SMU2ICEALBEDO_HPP */
