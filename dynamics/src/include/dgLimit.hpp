/*!
 * @file dgLimit.hpp
 * @date 30 April 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#ifndef __DGLIMIT_HPP
#define __DGLIMIT_HPP

#include "codeGenerationDGinGauss.hpp"
#include "dgVector.hpp"

namespace Nextsim {

//! Limit a dg vector from above
static void LimitMax(DGVector<1>& dg, double max) { dg.col(0) = dg.col(0).cwiseMin(max); }
static void LimitMax(DGVector<3>& dg, double max)
{
#pragma omp parallel for
    for (long int i = 0; i < dg.rows(); ++i) {
        dg(i, 0) = std::min(max, dg(i, 0));
        const double l0 = 2.0*std::max(fabs(dg(i, 1) + dg(i, 2)), fabs(dg(i, 1) - dg(i, 2)));
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

        double maxvalue  = (dg.block<1, 6>(i, 0) * PSI<6, 3>).maxCoeff();

	// part coming from the nonlinearities
        const double l0 = maxvalue - dg(i, 0);
	// value to big?
        if (maxvalue > max) {
	  // limit nonlinear part
            dg.block<1, 5>(i, 1) *= (max - dg(i, 0)) / l0;
        }
    }
}

//! Limit a dg vector from above in the gauss nodes
static void LimitMin(DGVector<1>& dg, double min) { dg.col(0) = dg.col(0).cwiseMax(min); }
static void LimitMin(DGVector<3>& dg, double min)
{
#pragma omp parallel for
    for (long int i = 0; i < dg.rows(); ++i) {
        dg(i, 0) = std::max(min, dg(i, 0));
        const double l0 = 2.0*std::max(fabs(dg(i, 1) + dg(i, 2)), fabs(dg(i, 1) - dg(i, 2)));
        if (l0 == 0)
            continue;
	
        const double ex = dg(i, 0) - l0 - min;
        if (ex < 0) {
            dg(i, 1) *= (dg(i, 0)-min) / l0;
            dg(i, 2) *= (dg(i, 0)-min) / l0;
        }
    }
}
static void LimitMin(DGVector<6>& dg, double min)
{
#pragma omp parallel for
    for (long int i = 0; i < dg.rows(); ++i) {
        dg(i, 0) = std::max(min, dg(i, 0));

        double minvalue = (dg.block<1, 6>(i, 0) * PSI<6, 3>).minCoeff();

	// part coming from nonlinearity
        const double l0 = minvalue - dg(i, 0);
        if (minvalue < min) {
            dg.block<1, 5>(i, 1) *= (dg(i, 0)-min) / l0;
        }
    }
}

} /* namespace Nextsim */

#endif /* __DGLIMIT_HPP */
