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

    template <int N>
    static DGVector<N>& hField2dg(const HField& h, DGVector<N>& dg)
    {
        dg.col(0) = h.data();
        return dg;
    }

    template <>
    inline DGVector<1>& hField2dg(const HField& h, DGVector<1>& dg)
    {
        return ma2dg(h, dg);
    }

    template <int N>
    static HField& dg2hField(const DGVector<N>& dg, HField& h)
    {
        h.setData(dg.col(0));
        return h;
    }

    template <>
    inline HField& dg2hField(const DGVector<1>& dg, HField& h)
    {
        return dg2ma(dg, h);
    }

};

} /* namespace Nextsim */

#endif /* DGMODELARRAY_HPP */
