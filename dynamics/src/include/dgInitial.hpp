/*!
 * @file dgInitial.hpp
 * @date July 10 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#ifndef __DGINITIAL_HPP
#define __DGINITIAL_HPP

#include "cgVector.hpp"
#include "dgVector.hpp"

#include "ParametricTools.hpp"

namespace Nextsim {

typedef std::function<double(double, double)> InitialOp;

struct SmoothInitial {
public:
    double operator()(double x, double y) const
    {
        double r2 = ((x - 0.7) * (x - 0.7) + (y - 0.7) * (y - 0.7));
        return exp(-200.0 * r2);
    }
};

struct PyramidInitial {
public:
    double operator()(double x, double y) const
    {
        double r2 = ((x - 0.3) * (x - 0.3) + (y - 0.6) * (y - 0.6));
        return std::max(0.0, 1.0 - 10.0 * sqrt(r2));
    }
};

struct BoxInitial {
public:
    double operator()(double x, double y) const
    {
        double r2 = ((x - 0.65) * (x - 0.65) + (y - 0.3) * (y - 0.3));
        if (sqrt(r2) < 0.1)
            return 1.0;
        return 0.0;
    }
};

struct MixedInitial {
public:
    double operator()(double x, double y) const
    {
        return SmoothInitial()(x, y) + BoxInitial()(x, y) + PyramidInitial()(x, y);
    }
};

////////////////////////////////////////////////// New Interface - ParametricMesh

////////////////////////////////////////////////// OLD Interfvace


} /* namespace Nextsim */

#endif /* __DGINITIAL_H */
