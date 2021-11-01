/*
 * @file IIceAlbedo.hpp
 *
 * @date Sep 28, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_IICEALBEDO_HPP
#define SRC_INCLUDE_IICEALBEDO_HPP

namespace Nextsim {
class IIceAlbedo {
public:
    virtual ~IIceAlbedo() = default;
    virtual double albedo(double temperature, double snowThickness) = 0;
};
}
#endif /* SRC_INCLUDE_IICEALBEDO_HPP */
