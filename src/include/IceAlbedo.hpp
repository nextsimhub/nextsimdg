/*
 * @file IceAlbedo.hpp
 *
 * @date Sep 22, 2021
 * @author Tim Spain
 */

#ifndef SRC_INCLUDE_ICEALBEDO_HPP
#define SRC_INCLUDE_ICEALBEDO_HPP

namespace Nextsim {

class IceAlbedo {
public:
    static double albedo(double temperature, double snowThickness);
};

}
#endif /* SRC_INCLUDE_ICEALBEDO_HPP */
