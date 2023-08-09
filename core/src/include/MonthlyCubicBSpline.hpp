/*!
 * @file MonthlyCubicBSpline.hpp
 *
 * @date Oct 14, 2022
 * @author Einar Örn Ólason <einar.olason@nersc.no>
 */

#ifndef MONTHLYCUBICBSPLINE_HPP
#define MONTHLYCUBICBSPLINE_HPP

#include <boost/version.hpp>

#if BOOST_VERSION >= 107200
#define NEW_SPLINES 1
#else
#define NEW_SPLINES 0
#endif

#if NEW_SPLINES
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
#else
#include <boost/math/interpolators/cubic_b_spline.hpp>
#endif

namespace Nextsim {
/*!
 * @brief A class for constructing cubib b-splines for equally spaced monthly data. The functions
 * ensure that the result is fully cyclic and takes leap years into account. The heavy lifting is
 * done by boost's cardinal_cubic_b_spline.
 */
class monthlyCubicBSpline {

public:
#if NEW_SPLINES
    typedef boost::math::interpolators::cardinal_cubic_b_spline<double> bSpline;
#else
    typedef boost::math::cubic_b_spline<double> bSpline;
#endif
    /*!
     * @brief The constructor for a monthlyCubicBSpline object.
     * @param f A vector of length 12 containing the monthly values for the spline
     */
    explicit monthlyCubicBSpline(const std::vector<double>& f)
    {
        // Create a cyclic B-spline using a fold with a "ghost border"
        /* We need k+1 ghost points, plus one for each of the derivatives at the end, where k is the
         * order of the spline
         * (here k=3). So, with STL's vector.insert() we need ghostWidth = k + 1 + 2 = 6 */
        const int ghostWidth = 6;
        std::vector<double> y = f;
        y.insert(y.begin(), f.end() - ghostWidth, f.end());
        y.insert(y.end(), f.begin(), f.begin() + ghostWidth);

        // Normalise the year from 0 to 1 (easier for leap years later on)
        const double h = 1. / 12.;
        const double t0 = h / 2. - ghostWidth * h;

        // Use boost!
        m_spline = std::make_shared<bSpline>(y.begin(), y.end(), t0, h);
    };

    /*!
     * @brief The operator for returning spline value at a given day of year
     * @param dayOfYear The day of year and fraction
     * @param isLeap Boolean to indicate if its a leap year
     * @return
     */
    double operator()(double dayOfYear, bool isLeap)
    {
        // Convert from day-of-year to normalised [0, 1]
        double fracYear;
        if (isLeap)
            fracYear = dayOfYear / 366.;
        else
            fracYear = dayOfYear / 365.;

        // Use boost!
        return (*m_spline)(fracYear);
    };

private:
    std::shared_ptr<bSpline> m_spline;
};
}

#endif // MONTHLYCUBICBSPLINE_HPP
