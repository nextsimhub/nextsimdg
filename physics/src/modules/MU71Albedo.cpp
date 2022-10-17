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

double MU71Albedo::albedo(double temperature, double snowThickness)
{
    // Fixed ice albedo
    if (snowThickness == 0.)
        return iceAlbedo;

    const double dayOfYear = M_tp.gmtime()->tm_yday;
    const bool isLeap = ((M_tp.gmtime()->tm_year % 4 == 0) && (M_tp.gmtime()->tm_year % 100 != 0))
        || (M_tp.gmtime()->tm_year % 400 == 0);

    return snowAlbedo(dayOfYear, isLeap);
}

}
