/*!
 * @file MU71Albedo.hpp
 *
 * @date Wed 24 Aug 2022 08:00:31 CEST
 * @author Einar Örn Ólason <einar.olason@nersc.no>
 */

#ifndef SEASONALICEALBEDO_HPP
#define SEASONALICEALBEDO_HPP

#include "IIceAlbedo.hpp"
#include "include/MonthlyCubicBSpline.hpp"
#include <vector>

namespace Nextsim {

/*!
 * @brief The implementation class for the a seasonal albedo following Maykut and Untersteiener's
 * (1971) table 1. Only useful for comparison with that paper and derived setups.
 */
class MU71Albedo : public IIceAlbedo {

public:
    MU71Albedo();

    /*!
     * @brief Returns the tabulated ice surface short wave albedo from Maykut and Untersteiner
     * (1971)
     * @param temperature The surface (ice or snow) temperature
     * @param snowThickness  Show thickness
     * @return
     */
    std::tuple<double, double> surfaceShortWaveBalance(
        double temperature, double snowThickness, double i0) override;

private:
    TimePoint M_tp;
    double dayOfYear;
    bool isLeap;

    void setTime(const TimePoint& tp) override
    {
        M_tp = tp;
        dayOfYear = M_tp.gmtime()->tm_yday;
        isLeap = ((M_tp.gmtime()->tm_year % 4 == 0) && (M_tp.gmtime()->tm_year % 100 != 0))
            || (M_tp.gmtime()->tm_year % 400 == 0);
    }

    // Monthly snow albedo from Maykut and Untersteiner (1971)
    const std::vector<double> albedoTable
        = { 0.85, 0.85, 0.83, 0.81, 0.82, 0.78, 0.64, 0.69, 0.84, 0.85, 0.85, 0.85 };
    //      Jan,  Feb,  Mar,  Apr,  Mai,  Jun,  Jul,  Aug,  Sept, Oct,  Nov,  Dec

    monthlyCubicBSpline snowAlbedo;
    const double iceAlbedo = 0.64;
};

}

#endif /* SEASONALICEALBEDO_HPP */
