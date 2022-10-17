//
// Created by Einar Ólason on 14/10/2022.
//

#ifndef MONTHLYCUBICBSPLINE_HPP
#define MONTHLYCUBICBSPLINE_HPP

#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>

namespace Nextsim {
class monthlyCubicBSpline {

public:
    explicit monthlyCubicBSpline(const std::vector<double>& f)
    {
        // Create a cyclic B-spline using a fold with a "ghost border"
        const int ghostWidth = 5;
        auto y = f;
        y.insert(y.begin(), f.end() - ghostWidth, f.end());
        y.insert(y.end(), f.begin(), f.begin() + ghostWidth);

        // Normalise the year from 0 to 1 (easier for leap years later on)
        const double h = 1. / 12.;
        const double t0 = h / 2. - ghostWidth * h;

        // Estimate the derivatives at the end of the "ghost border"
        const double start_deriv = (*(f.end() - ghostWidth) - *(f.end() - ghostWidth - 1)) / h;
        const double end_deriv = (*(f.begin() + ghostWidth) - *(f.begin() + ghostWidth - 1)) / h;

        // Use boost!
        m_spline = std::make_shared<boost::math::interpolators::cardinal_cubic_b_spline<double>>(
            y.begin(), y.end(), t0, h, start_deriv, end_deriv);
    };

    double operator()(double dayOfYear, bool isLeap = false)
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
    std::shared_ptr<boost::math::interpolators::cardinal_cubic_b_spline<double>> m_spline;
};
}

#endif // NEXTSIM_DG_MONTHLYCUBICBSPLINE_HPP
