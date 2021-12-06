/*!
 * @file SMU2IceAlbedo.hpp
 *
 * @date Sep 28, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_SMU2ICEALBEDO_HPP
#define SRC_INCLUDE_SMU2ICEALBEDO_HPP

#include "IIceAlbedo.hpp"

namespace Nextsim {

class SMU2IceAlbedo : public IIceAlbedo {
    double albedo(double temperature, double snowThickness);
};

}

#endif /* SRC_INCLUDE_SMU2ICEALBEDO_HPP */
