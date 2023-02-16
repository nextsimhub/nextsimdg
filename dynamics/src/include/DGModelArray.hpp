/*!
 * @file DGModelArray.hpp
 *
 * @date Oct 6, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef DGMODELARRAY_HPP
#define DGMODELARRAY_HPP

#include "include/ModelArray.hpp"
#include "dgVector.hpp"

#include <cassert>

namespace Nextsim {

class DGModelArray {
public:
    template <int N>
    static DGVector<N>& ma2dg(const ModelArray& ma, DGVector<N>& dg)
    {
        if (N == ma.components(0).size()) {
//            assert(N == ma.components(0).size());
            dg = ma.data().matrix();
        } else {
            // Assign only to the 0 component.
            dg(Eigen::all, 0) = ma.data().matrix();
        }
        return dg;
    }

    template <int N>
    static ModelArray& dg2ma(const  DGVector<N>& dg, ModelArray& ma)
    {
        if (N == ma.components(0).size()) {
            ma.setData(dg.data());
        } else {
            /* Assign the zero component as data. Since the setData function
             * takes a pointer to continuous data, the data needs to be copied
             * from the DGVector initially.
             */
            Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::ColMajor> dg0Data = dg(Eigen::all, 0);
            ma.setData(dg0Data.data());
        }
        return ma;
    }

};

} /* namespace Nextsim */

#endif /* DGMODELARRAY_HPP */
