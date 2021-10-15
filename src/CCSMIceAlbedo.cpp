/*
 * @file CCSMIceAlbedo.cpp
 *
 * @date Sep 22, 2021
 * @author Tim Spain
 */

#include "include/CCSMIceAlbedo.hpp"

#include <cmath>

/* Albedo scheme from ccsm 3 */
/* The scheme is simplified using the assumption that visible solar
 * radiation accounts for 52% of the spectrum and the near-infrared for
 * 48% (the same split as used in cice when running in stand-alone
 * mode). */

/* This is the ccsm3 scheme when alb_ice = 0.538 and alb_sn = 0.8256 */

namespace Nextsim {

const double ICE_ALBEDO = 0.538;
const double SNOW_ALBEDO = 0.8256;

double CCSMIceAlbedo::albedo(double temperature, double snowThickness)
{
    const double tLimit = -1.;
    double iceAlbedo = ICE_ALBEDO - std::fmax(0., 0.075 * (temperature - tLimit));
    double snowAlbedo = SNOW_ALBEDO - std::fmax(0., 0.124 * (temperature - tLimit));

    double snowCoverFraction = snowThickness/(snowThickness + 0.02);

    return snowCoverFraction * snowAlbedo + (1 - snowCoverFraction) * iceAlbedo;
}
}
