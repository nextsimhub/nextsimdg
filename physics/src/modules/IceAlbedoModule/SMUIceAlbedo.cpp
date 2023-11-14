/*!
 * @file SMUIceAlbedo.cpp
 *
 * @date Sep 22, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/SMUIceAlbedo.hpp"

namespace Nextsim {

/* This scheme mimics Semtner 76 and Maykut and Untersteiner 71 when
 * alb_ice = 0.64 and alb_sn = 0.85 */

const double ICE_ALBEDO = 0.64;
const double SNOW_ALBEDO = 0.85;

std::tuple<double, double> SMUIceAlbedo::albedo(double temperature, double snowThickness, double i0)
{
    double albedo, penSW;
    if (snowThickness > 0.) {
        albedo = SNOW_ALBEDO;
        penSW = 0.;
    } else {
        albedo = ICE_ALBEDO;
        penSW = i0;
    }
    return {albedo, penSW};
}
}
