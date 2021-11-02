/*
 * @file SMUIceAlbedo.hpp
 *
 * @date Sep 28, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_SMUICEALBEDO_HPP
#define SRC_INCLUDE_SMUICEALBEDO_HPP

#include "IIceAlbedo.hpp"

namespace Nextsim {

class SMUIceAlbedo : public IIceAlbedo {
    double albedo(double temperature, double snowThickness);
};

}

#endif /* SRC_INCLUDE_SMUICEALBEDO_HPP */
