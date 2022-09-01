/*!
 * @file MonthlyIceAlbedo.hpp
 *
 * @date Wed 24 Aug 2022 08:00:31 CEST
 * @author Einar Örn Ólason <einar.olason@nersc.no>
 */

#ifndef SRC_INCLUDE_SEASONALICEALBEDO_HPP
#define SRC_INCLUDE_SEASONALICEALBEDO_HPP

#include "include/TableLookup.hpp"
#include "IIceAlbedo.hpp"
#include <vector>

namespace Nextsim {

//! The implementation class for the a seasonal albedo following Maykut and Untersteiener's (1971)
//! table 1. Only useful for comparison with that paper and derived setups.
class MonthlyIceAlbedo : public IIceAlbedo {

public:
    /*!
     * @brief Returns the tabulated ice surface short wave albedo from Maykut and Untersteiner
     * (1971)
     */
    double albedo(double temperature, double snowThickness);

private:
    TimePoint M_tp;
    void setTime(const TimePoint& tp) { M_tp = tp;}
};

}

#endif /* SRC_INCLUDE_SEASONALICEALBEDO_HPP */
