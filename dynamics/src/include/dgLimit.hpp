/*!
 * @file dgLimit.hpp
 * @date 30 April 2022
 * @author Thomas Richter <thomas.richter@ovgu.no>
 */

#ifndef __DGLIMIT_HPP
#define __DGLIMIT_HPP

#include "codeGenerationDGinGauss.hpp"
#include "dgVector.hpp"

namespace Nextsim {

//! Limit a dg vector from above
static void LimitMax(DGVector<1>& dg, double max)
{
    dg.col(0) = dg.col(0).cwiseMin(max);
}
static void LimitMax(DGVector<3>& dg, double max)
{
#pragma omp parallel for
    for (long int i = 0; i < dg.rows(); ++i) {
        dg(i, 0) = std::min(max, dg(i, 0));
        const double l0 = std::max(fabs(dg(i, 1) + dg(i, 2)), fabs(dg(i, 1) - dg(i, 2)));
        if (l0 == 0)
            continue;
        const double ex = dg(i, 0) + l0 - max;
        if (ex > 0) {
            dg(i, 1) *= (max - dg(i, 0)) / l0;
            dg(i, 2) *= (max - dg(i, 0)) / l0;
        }
    }
}

//! Limit a dg vector from above in the gauss nodes
static void LimitMax(DGVector<6>& dg, double max)
{
#pragma omp parallel for
    for (long int i = 0; i < dg.rows(); ++i) {
        dg(i, 0) = std::min(max, dg(i, 0));

        LocalDGVector<9> dgingauss = dg.block<1, 6>(i, 0) * PSI<6, 3>;
        const double l0 = dgingauss.maxCoeff() - dg(i, 0);
        if (l0 > max - dg(i, 0)) {
            dg.block<1, 5>(i, 1) *= (max - dg(i, 0)) / l0;
        }
    }
}


//! Limit a dg vector from above in the gauss nodes
static void LimitMin(DGVector<1>& dg, double min)
{
    dg.col(0) = dg.col(0).cwiseMax(min);
}
static void LimitMin(DGVector<3>& dg, double min)
{
#pragma omp parallel for
    for (long int i = 0; i < dg.rows(); ++i) {
        dg(i, 0) = std::max(min, dg(i, 0));
        const double l0 = std::max(fabs(dg(i, 1) + dg(i, 2)), fabs(dg(i, 1) - dg(i, 2)));
        if (l0 == 0)
            continue;
        const double ex = dg(i, 0) - l0 - min;
        if (ex < 0) {
            dg(i, 1) *= (min - dg(i, 0)) / l0;
            dg(i, 2) *= (min - dg(i, 0)) / l0;
        }
    }
}
static void LimitMin(DGVector<6>& dg, double min)
{
#pragma omp parallel for
    for (long int i = 0; i < dg.rows(); ++i) {
        dg(i, 0) = std::max(min, dg(i, 0));

        LocalDGVector<9> dgingauss = dg.block<1, 6>(i, 0) * PSI<6, 3>;
        const double l0 = dgingauss.minCoeff() - dg(i, 0);
        if (l0 < min - dg(i, 0)) {
            dg.block<1, 5>(i, 1) *= (min - dg(i, 0)) / l0;
        }
    }
}

} /* namespace Nextsim */

#endif /* __DGLIMIT_HPP */
