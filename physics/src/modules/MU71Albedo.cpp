/*!
 * @file MU71Albedo.cpp
 *
 * @date Wed 24 Aug 2022 08:05:49 CEST
 * @author Einar Örn Ólason <einar.olason@nersc.no>
 */

#include "include/MU71Albedo.hpp"

namespace Nextsim {
double MU71Albedo::albedo(double temperature, double snowThickness)
{
    // Monthly albedos from Maykut and Untersteiner (1971)
    const std::vector<double> albedoTable
        = { 0.85, 0.85, 0.83, 0.81, 0.82, 0.78, 0.64, 0.69, 0.84, 0.85, 0.85, 0.85 };
    //      Jan,  Feb,  Mar,  Apr,  Mai,  Jun,  Jul,  Aug,  Sept, Oct,  Nov,  Dec

    const double dayOfYear = M_tp.gmtime()->tm_yday;
    const bool isLeap = ((M_tp.gmtime()->tm_year % 4 == 0) && (M_tp.gmtime()->tm_year % 100 != 0))
        || (M_tp.gmtime()->tm_year % 400 == 0);

    return TableLookup::monthlyLinearLUT(albedoTable, dayOfYear, isLeap);
}

}
