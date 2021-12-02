/*!
 * @file CCSMIceAlbedo.hpp
 *
 * @date Sep 28, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_CCSMICEALBEDO_HPP
#define SRC_INCLUDE_CCSMICEALBEDO_HPP

#include "Configured.hpp"
#include "IIceAlbedo.hpp"

namespace Nextsim {

class CCSMIceAlbedo : public IIceAlbedo, public Configured<CCSMIceAlbedo> {
public:
    double albedo(double temperature, double snowThickness) override;

    void configure() override;
private:
    static double iceAlbedo;
    static double snowAlbedo;
};

}

#endif /* SRC_INCLUDE_CCSMICEALBEDO_HPP */
