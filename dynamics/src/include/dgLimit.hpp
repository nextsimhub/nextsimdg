/*!
 * @author  Thomas Richter <thomas.richter@ovgu.de>
 */

#ifndef __DGLIMIT_HPP
#define __DGLIMIT_HPP

#include "codeGenerationDGinGauss.hpp"
#include "dgVector.hpp"

namespace Nextsim {

/*!
 * Computes the minimum and maximum of the dQ1 polynomial
 *   f(X,Y) = c0 + c1*X + c2*Y + c3*X*Y
 * over the reference cell X,Y in [-1/2,1/2]. f is linear along each edge,
 * so the extrema are attained in the corners.
 */
static void DG4Extrema(const DGVector<4>& dg, long int i, FloatType& minval, FloatType& maxval)
{
    const FloatType v0 = dg(i, 0) - 0.5_ft * dg(i, 1) - 0.5_ft * dg(i, 2) + 0.25_ft * dg(i, 3);
    const FloatType v1 = dg(i, 0) + 0.5_ft * dg(i, 1) - 0.5_ft * dg(i, 2) - 0.25_ft * dg(i, 3);
    const FloatType v2 = dg(i, 0) - 0.5_ft * dg(i, 1) + 0.5_ft * dg(i, 2) - 0.25_ft * dg(i, 3);
    const FloatType v3 = dg(i, 0) + 0.5_ft * dg(i, 1) + 0.5_ft * dg(i, 2) + 0.25_ft * dg(i, 3);
    minval = std::min(std::min(v0, v1), std::min(v2, v3));
    maxval = std::max(std::max(v0, v1), std::max(v2, v3));
}

/*!
 * Computes the minimum and maximum of the dQ2 polynomial
 *   f(X,Y) = c0 + c1*X + c2*Y + c3*X*Y + c4*qX + c5*qY + c6*qX*Y + c7*X*qY + c8*qX*qY
 * with qX = X^2-1/12, qY = Y^2-1/12, over the reference cell X,Y in [-1/2,1/2].
 * On each line Y = const (resp. X = const) f is a 1D quadratic, whose extremum on
 * [-1/2,1/2] is taken exactly (vertex of the parabola or end points).
 * 1. NLINES+1 equidistant lines in each direction, the edges included: exact on the edges
 *    and in the corners.
 * 2. The interior critical point is then refined from the best point of 1. by alternating
 *    exact line searches in X and Y (the value can only improve).
 */
static void DG9Extrema(const DGVector<9>& dg, long int i, FloatType& minval, FloatType& maxval)
{
    const Eigen::Matrix<FloatType, 1, 9> c = dg.row(i);
    auto f = [&](FloatType X, FloatType Y) {
        const FloatType qX = X * X - 1.0_ft / 12.0_ft;
        const FloatType qY = Y * Y - 1.0_ft / 12.0_ft;
        return c(0) + c(1) * X + c(2) * Y + c(3) * X * Y + c(4) * qX + c(5) * qY + c(6) * qX * Y
            + c(7) * X * qY + c(8) * qX * qY;
    };

    // best X on the line Y = s (largest value if sign = 1, smallest if sign = -1)
    auto bestX = [&](FloatType s, FloatType sign) {
        const FloatType qs = s * s - 1.0_ft / 12.0_ft;
        const FloatType B = c(1) + c(3) * s + c(7) * qs; // f = A + B X + C X^2
        const FloatType C = c(4) + c(6) * s + c(8) * qs;
        FloatType X = (sign * f(0.5_ft, s) > sign * f(-0.5_ft, s)) ? 0.5_ft : -0.5_ft;
        if (C != 0.0) {
            const FloatType Xc = -B / (2.0_ft * C);
            if (Xc >= -0.5_ft && Xc <= 0.5_ft && sign * f(Xc, s) > sign * f(X, s))
                X = Xc;
        }
        return X;
    };
    // best Y on the line X = s
    auto bestY = [&](FloatType s, FloatType sign) {
        const FloatType qs = s * s - 1.0_ft / 12.0_ft;
        const FloatType B = c(2) + c(3) * s + c(6) * qs; // f = A + B Y + C Y^2
        const FloatType C = c(5) + c(7) * s + c(8) * qs;
        FloatType Y = (sign * f(s, 0.5_ft) > sign * f(s, -0.5_ft)) ? 0.5_ft : -0.5_ft;
        if (C != 0.0) {
            const FloatType Yc = -B / (2.0_ft * C);
            if (Yc >= -0.5_ft && Yc <= 0.5_ft && sign * f(s, Yc) > sign * f(s, Y))
                Y = Yc;
        }
        return Y;
    };

    // sign = 1: maximum, sign = -1: minimum
    auto extremum = [&](FloatType sign) {
        FloatType Xb = -0.5_ft, Yb = -0.5_ft;
        FloatType vb = sign * f(Xb, Yb);
        auto consider = [&](FloatType X, FloatType Y) {
            const FloatType v = sign * f(X, Y);
            if (v > vb) {
                vb = v;
                Xb = X;
                Yb = Y;
            }
        };

        constexpr int NLINES = 8;
        for (int k = 0; k <= NLINES; ++k) {
            const FloatType s = -0.5_ft + FloatType(k) / NLINES;
            consider(bestX(s, sign), s); // line Y = s
            consider(s, bestY(s, sign)); // line X = s
        }

        constexpr int NREFINE = 50;
        for (int it = 0; it < NREFINE; ++it) {
            consider(bestX(Yb, sign), Yb);
            consider(Xb, bestY(Xb, sign));
        }
        return sign * vb;
    };

    maxval = extremum(1.0_ft);
    minval = extremum(-1.0_ft);
}

//! Limit a dg vector from above
static void LimitMax(DGVector<1>& dg, FloatType max) { dg.col(0) = dg.col(0).cwiseMin(max); }
static void LimitMax(DGVector<3>& dg, FloatType max)
{
#pragma omp parallel for
    for (long int i = 0; i < dg.rows(); ++i) {
        dg(i, 0) = std::min(max, dg(i, 0));
        const FloatType l0
            = 2.0_ft * std::max(fabs(dg(i, 1) + dg(i, 2)), fabs(dg(i, 1) - dg(i, 2)));
        if (l0 == 0)
            continue;
        const FloatType ex = dg(i, 0) + l0 - max;
        if (ex > 0) {
            dg(i, 1) *= (max - dg(i, 0)) / l0;
            dg(i, 2) *= (max - dg(i, 0)) / l0;
        }
    }
}

//! Limit a dg vector from above in the corners (exact for the bilinear polynomial)
static void LimitMax(DGVector<4>& dg, FloatType max)
{
#pragma omp parallel for
    for (long int i = 0; i < dg.rows(); ++i) {
        dg(i, 0) = std::min(max, dg(i, 0));

        FloatType minvalue, maxvalue;
        DG4Extrema(dg, i, minvalue, maxvalue);
        // part coming from the nonlinearities
        const FloatType l0 = maxvalue - dg(i, 0);
        if (maxvalue > max) {
            dg.block<1, 3>(i, 1) *= (max - dg(i, 0)) / l0;
        }
    }
}

//! Limit a dg vector from above in the gauss nodes
static void LimitMax(DGVector<6>& dg, FloatType max)
{
#pragma omp parallel for
    for (long int i = 0; i < dg.rows(); ++i) {
        dg(i, 0) = std::min(max, dg(i, 0));

        FloatType maxvalue = (dg.block<1, 6>(i, 0) * PSI<6, 3>).maxCoeff();

        // part coming from the nonlinearities
        const FloatType l0 = maxvalue - dg(i, 0);
        // value to big?
        if (maxvalue > max) {
            // limit nonlinear part
            dg.block<1, 5>(i, 1) *= (max - dg(i, 0)) / l0;
        }
    }
}

//! Limit a dg vector from above with the extrema over the cell (DG9Extrema)
static void LimitMax(DGVector<9>& dg, FloatType max)
{
#pragma omp parallel for
    for (long int i = 0; i < dg.rows(); ++i) {
        dg(i, 0) = std::min(max, dg(i, 0));

        FloatType minvalue, maxvalue;
        DG9Extrema(dg, i, minvalue, maxvalue);
        // part coming from the nonlinearities
        const FloatType l0 = maxvalue - dg(i, 0);
        if (maxvalue > max) {
            dg.block<1, 8>(i, 1) *= (max - dg(i, 0)) / l0;
        }
    }
}

//! Limit a dg vector from above in the gauss nodes
static void LimitMin(DGVector<1>& dg, FloatType min) { dg.col(0) = dg.col(0).cwiseMax(min); }
static void LimitMin(DGVector<3>& dg, FloatType min)
{
#pragma omp parallel for
    for (long int i = 0; i < dg.rows(); ++i) {
        dg(i, 0) = std::max(min, dg(i, 0));
        const FloatType l0
            = 2.0_ft * std::max(fabs(dg(i, 1) + dg(i, 2)), fabs(dg(i, 1) - dg(i, 2)));
        if (l0 == 0)
            continue;

        const FloatType ex = dg(i, 0) - l0 - min;
        if (ex < 0) {
            dg(i, 1) *= (dg(i, 0) - min) / l0;
            dg(i, 2) *= (dg(i, 0) - min) / l0;
        }
    }
}
static void LimitMin(DGVector<4>& dg, FloatType min)
{
#pragma omp parallel for
    for (long int i = 0; i < dg.rows(); ++i) {
        dg(i, 0) = std::max(min, dg(i, 0));

        FloatType minvalue, maxvalue;
        DG4Extrema(dg, i, minvalue, maxvalue);
        // part coming from the nonlinearities
        const FloatType l0 = minvalue - dg(i, 0);
        if (minvalue < min) {
            dg.block<1, 3>(i, 1) *= (min - dg(i, 0)) / l0;
        }
    }
}
static void LimitMin(DGVector<6>& dg, FloatType min)
{
#pragma omp parallel for
    for (long int i = 0; i < dg.rows(); ++i) {
        dg(i, 0) = std::max(min, dg(i, 0));

        FloatType minvalue = (dg.block<1, 6>(i, 0) * PSI<6, 3>).minCoeff();

        // part coming from nonlinearity
        const FloatType l0 = minvalue - dg(i, 0);
        if (minvalue < min) {
            dg.block<1, 5>(i, 1) *= (min - dg(i, 0)) / l0;
        }
    }
}
static void LimitMin(DGVector<9>& dg, FloatType min)
{
#pragma omp parallel for
    for (long int i = 0; i < dg.rows(); ++i) {
        dg(i, 0) = std::max(min, dg(i, 0));

        FloatType minvalue, maxvalue;
        DG9Extrema(dg, i, minvalue, maxvalue);
        // part coming from the nonlinearities
        const FloatType l0 = minvalue - dg(i, 0);
        if (minvalue < min) {
            dg.block<1, 8>(i, 1) *= (min - dg(i, 0)) / l0;
        }
    }
}

} /* namespace Nextsim */

#endif /* __DGLIMIT_HPP */
