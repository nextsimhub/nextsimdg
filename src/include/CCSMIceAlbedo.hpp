/*
 * @file CCSMIceAlbedo.hpp
 *
 * @date Sep 28, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_CCSMICEALBEDO_HPP
#define SRC_INCLUDE_CCSMICEALBEDO_HPP

#include "IICeAlbedo.hpp"

namespace Nextsim {

class CCSMIceAlbedo : public IIceAlbedo {
    double albedo(double temperature, double snowThickness);
};

}

#endif /* SRC_INCLUDE_CCSMICEALBEDO_HPP */
