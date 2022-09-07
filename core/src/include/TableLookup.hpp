//
// Created by Einar Ólason on 29/08/2022.
//

#ifndef NEXTSIM_DG_TABLELOOKUP_HPP
#define NEXTSIM_DG_TABLELOOKUP_HPP

#include <vector>

namespace Nextsim {

/*!
 * A small class containing static lookup-table functions.
 */
class TableLookup {

public:
    // A general table lookup with a linear interpolation between tabulated values
    /*!
     * @brief A general table lookup with a linear interpolation between tabulated values
     * @param x: Input x-values
     * @param y: Input y-values
     * @param xx: The value on x that y should be interpolated to
     * @param cyclic: Boolean option to treat x and y as cyclic values
     * @return The interpolated y value at xx
     */
    double static linearLUT(const std::vector<double>& x, const std::vector<double>& y,
        const double xx, const bool cyclic = false)
    {
        // The standard case for lower and upper bound
        size_t i = std::lower_bound(x.begin(), x.end(), xx) - x.begin();
        size_t im1 = i - 1;

        // The cyclic cases
        if (cyclic) {
            if (i == x.size())
                i = 0;
            else if (i == 0)
                im1 = x.size() - 1;
        }

        return y[im1] + (xx - x[im1]) * (y[i] - y[im1]) / (x[i] - x[im1]);
    }

    // A general table lookup for mid-monthly values
public:
    /*!
     * @brief A general table lookup with a linear interpolation between tabulated monthly values
     * @param y: Input tabulated values
     * @param dayOfYear: The day on which the tabulated value should be interpolated
     * @return The interpolated value at dayOfYear
     */
    double static monthlyLinearLUT(
        const std::vector<double>& y, double dayOfYear, bool isLeap = false)
    {
        //                                1   2   3   4   5   6   7   8   9  10  11  12
        std::vector<int> monthLength = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
        if (isLeap)
            ++monthLength[1];

        // We use the mid-month day-of-year as a double for the interpolation
        std::vector<double> midMonth(monthLength.size());
        double monthEnd = 0;
        for (int i = 0; i < monthLength.size(); ++i) {
            midMonth[i] = double(monthEnd) + double(monthLength[i]) / 2.;
            monthEnd += monthLength[i];
        }

        return linearLUT(midMonth, y, dayOfYear, true);
    }
};

}

#endif // NEXTSIM_DG_TABLELOOKUP_HPP
