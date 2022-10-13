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
        assert(N == ma.components(0).size());
        return dg = ma.data().matrix();
    }

    template <int N>
    static ModelArray& dg2ma(const  DGVector<N>& dg, ModelArray& ma)
    {
        assert(N == ma.components(0).size());
        ma.setData(dg.data());
        return ma;
    }

};

} /* namespace Nextsim */

#endif /* DGMODELARRAY_HPP */
