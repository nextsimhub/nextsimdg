/*!
 * @file CGModelArray.hpp
 *
 * @date Oct 13, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CGMODELARRAY_HPP
#define CGMODELARRAY_HPP

#include "cgVector.hpp"
#include "include/ModelArray.hpp"

namespace Nextsim {

class CGModelArray {
public:
    template <int CG> static CGVector<CG>& ma2cg(const ModelArray& ma, CGVector<CG>& cg)
    {
        cg = ma.data().matrix();
        return cg;
    }

    template <int CG> static ModelArray& cg2ma(const CGVector<CG>& cg, ModelArray& ma)
    {
        ma.setData(cg.data());
        return ma;
    }

    template <int CG> static ModelArray::MultiDim cgDimensions(const ModelArray::MultiDim& hDims)
    {
        ModelArray::MultiDim cgDims(hDims);
        for (size_t i = 0; i < cgDims.size(); ++i) {
            cgDims[i] *= CG;
            ++cgDims[i];
        }
        return cgDims;
    }
};
}

#endif /* CGMODELARRAY_HPP */
