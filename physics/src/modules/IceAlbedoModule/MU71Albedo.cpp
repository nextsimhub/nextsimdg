/*!
 * @file MU71Albedo.cpp
 *
 * @date Wed 24 Aug 2022 08:05:49 CEST
 * @author Einar Örn Ólason <einar.olason@nersc.no>
 */

#include "include/MU71Albedo.hpp"

namespace Nextsim {

MU71Albedo::MU71Albedo()
    : snowAlbedo(albedoTable)
{
}

std::tuple<double, double> MU71Albedo::surfaceShortWaveBalance(
    double temperature, double snowThickness, double i0)
{
    double albedo, penSW;

    // Fixed ice albedo
    if (snowThickness == 0.) {
        albedo = iceAlbedo;
        penSW = i0;
    } else {
        albedo = snowAlbedo(dayOfYear, isLeap);
        penSW = 0.;
    }

    return std::make_tuple(albedo, penSW);
}

}
