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

double SMUIceAlbedo::albedo(double temperature, double snowThickness)
{
    if (snowThickness > 0.) {
        return SNOW_ALBEDO;
    } else {
        return ICE_ALBEDO + 0.4 * (1 - ICE_ALBEDO) * 0.17;// FIXME NextsimPhysics::i0();
    }
}
}
