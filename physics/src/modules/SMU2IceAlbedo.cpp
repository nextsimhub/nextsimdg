/*!
 * @file SMU2IceAlbedo.cpp
 *
 * @date Sep 22, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/SMU2IceAlbedo.hpp"
#include "include/NextsimPhysics.hpp"

#include <cmath>

namespace Nextsim {

/* This scheme mimics Semtner 76 and Maykut and Untersteiner 71 when
 * alb_ice = 0.64 and alb_sn = 0.85 */

const double ICE_ALBEDO = 0.64;
const double SNOW_ALBEDO = 0.85;

double SMU2IceAlbedo::albedo(double temperature, double snowThickness)
{
    if (snowThickness > 0.) {
        return std::fmin(
            SNOW_ALBEDO, ICE_ALBEDO + (SNOW_ALBEDO - ICE_ALBEDO) * snowThickness / 0.2);
    } else {
        return ICE_ALBEDO + 0.4 * (1 - ICE_ALBEDO) * 0.17;// FIXME NextsimPhysics::i0();
    }
}
}
