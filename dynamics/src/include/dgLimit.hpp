/*!
 * @file MeshdgLimit.hpp
 * @date 30 April 2022
 * @author Thomas Richter <thomas.richter@ovgu.no>
 */

#ifndef __DGLIMIT_HPP
#define __DGLIMIT_HPP

#include "codeGenerationDGinGauss.hpp"
#include "dgVector.hpp"

namespace Nextsim {

//! Limit a dg vector from above
void LimitMax(CellVector<1>& dg, double max)
{
    dg.col(0) = dg.col(0).cwiseMin(max);
}
void LimitMax(CellVector<3>& dg, double max)
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

// limits in the gauss nodes
void LimitMax(CellVector<6>& dg, double max)
{
#pragma omp parallel for
    for (long int i = 0; i < dg.rows(); ++i) {
        dg(i, 0) = std::min(max, dg(i, 0));

        LocalCellVector<9> dgingauss = dg.block<1, 6>(i, 0) * PSI<6,3>;
        const double l0 = dgingauss.maxCoeff() - dg(i, 0);
        if (l0 > max - dg(i, 0)) {
            dg.block<1, 5>(i, 1) *= (max - dg(i, 0)) / l0;
        }
    }
}

//   //! Limit a dg vector from above
//   void LimitMin(CellVector<1>& dg, double min)
//   {
//     for (long int i = 0; i < dg.rows(); ++i) {
//         dg(i, 0) = std::min(max, dg(i, 0));
//         const double l0 = std::max(fabs(dg(i, 1) + dg(i, 2)), fabs(dg(i, 1) - dg(i, 2)));
//         if (l0 == 0)
//             continue;
//         const double ex = dg(i, 0) + l0 - max;
//         if (ex > 0) {
//             dg(i, 1) *= (max - dg(i, 0)) / l0;
//             dg(i, 2) *= (max - dg(i, 0)) / l0;
//         }
//     }
// }

//! Limit a dg vector from above
void LimitMin(CellVector<1>& dg, double min)
{
    dg.col(0) = dg.col(0).cwiseMax(min);
}
void LimitMin(CellVector<3>& dg, double min)
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
void LimitMin(CellVector<6>& dg, double min)
{
#pragma omp parallel for
    for (long int i = 0; i < dg.rows(); ++i) {
        dg(i, 0) = std::max(min, dg(i, 0));

        LocalCellVector<9> dgingauss = dg.block<1, 6>(i, 0) * PSI<6,3>;
        const double l0 = dgingauss.minCoeff() - dg(i, 0);
        if (l0 < min - dg(i, 0)) {
            dg.block<1, 5>(i, 1) *= (min - dg(i, 0)) / l0;
        }
    }
}

//     for (long int i = 0; i < dg.rows(); ++i) {
//         dg(i, 0) = std::max(min, dg(i, 0));
//         const double l0 = std::max(fabs(dg(i, 1) + dg(i, 2)), fabs(dg(i, 1) - dg(i, 2)));
//         if (l0 == 0)
//             continue;
//         const double ex = dg(i, 0) - l0 - min;
//         if (ex < 0) {
//             dg(i, 1) *= (min - dg(i, 0)) / l0;
//             dg(i, 2) *= (min - dg(i, 0)) / l0;
//         }
//     }
// }
// >>>>>>> e77045893a2f9b81c3242f0d4398b00f3b16af1b

} /* namespace Nextsim */

#endif /* __DGLIMIT_HPP */
