/*!
 * @file DGModelArray.hpp
 *
 * @date Oct 6, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef DGMODELARRAY_HPP
#define DGMODELARRAY_HPP

#include "dgVector.hpp"
#include "include/ModelArray.hpp"

#include <cassert>

namespace Nextsim {

class DGModelArray {
public:
    template <int N> static DGVector<N>& ma2dg(const ModelArray& ma, DGVector<N>& dg)
    {
        if (N == ma.components(0).size()) {
            //            assert(N == ma.components(0).size());
            dg = ma.data().matrix();
        } else {
            // Assign only to the 0 component.
            // dg(Eigen::all, 0) = ma.data().matrix();
            dg.col(0) = ma.data().matrix();
        }
        return dg;
    }

    /*!
     * Copies the contents of a 2D slice from a ModelArray to a DGVector, where the Model Array
     * has more than 2 dimensions.
     *
     * @param ma The source ModelArray, with more than 2 dimensions
     * @param dg The target DGVector
     * @param kField The selected z-level of ModelArray ma.
     */
    template <int DG> static DGVector<DG>& ma2dg(const ModelArray& ma, DGVector<DG>& dg, size_t kField)
    {
        // Use the 2D optimized version if possible
        if (ma.nDimensions() == 2) return ma2dg(ma, dg);

        // Extract the desired level in the first dimension
        // The size of a 2D domain is taken to be the product of the first two dimensions of ma
        ModelArray::MultiDim maDims = ma.dimensions();
        size_t size2D = maDims[0] * maDims[1];
        // With DG components
        if (DG == ma.components(0).size()) {
            dg = ma.data().matrix()(Eigen::seqN(kField * size2D, size2D), Eigen::all);
        } else {
            dg.col(0) = ma.data().matrix()(Eigen::seqN(kField * size2D, size2D), 0);
        }
        return dg;
    }

    /*!
     * Copies the contents of 2D slice from a ModelArray to a DGVector, where the Model Array
     * has more than 2 dimensions.
     *
     * @param ma The source ModelArray, with more than 2 dimensions
     * @param dg The target DGVector
     * @param sliceIndices The selected z-level of ModelArray ma.
     */
    template <int DG> static DGVector<DG>& ma2dg(const ModelArray& ma, DGVector<DG>& dg, const ModelArray::MultiDim& sliceIndices)
    {
        // Use the 2D optimized version if possible
        if (ma.nDimensions() == 2) return ma2dg(ma, dg);

        return ma2dg(ma, dg, kFromLocation(sliceIndices, ma));
    }

    template <int N> static ModelArray& dg2ma(const DGVector<N>& dg, ModelArray& ma)
    {
        if (N == ma.components(0).size()) {
            ma.setData(dg.data());
        } else {
            /* Assign the zero component as data. Since the setData function
             * takes a pointer to continuous data, the data needs to be copied
             * from the DGVector initially.
             */
            // Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::ColMajor> dg0Data = dg(Eigen::all,
            // 0); ma.setData(dg0Data.data());
            ma.setData(dg.col(0));
        }
        return ma;
    }

    template <int DG> static ModelArray& dg2ma(const DGVector<DG>& dg, ModelArray& ma, size_t kIndex)
    {
        // Use the 2D optimized version if possible
        if (ma.nDimensions() == 2) return dg2ma(dg, ma);

        // Extract the desired level in the first dimension
        // The size of a 2D domain is taken to be the product of the first two dimensions of ma
        ModelArray::MultiDim maDims = ma.dimensions();
        size_t size2D = maDims[0] * maDims[1];
        if (DG == ma.components(0).size()) {
            ma.setData(dg.array(), kIndex* size2D, size2D);
        } else {
            ma.setData(dg.array().col(0), kIndex* size2D, size2D);
        }
        return ma;
    }

    template <int DG> static ModelArray& dg2ma(const DGVector<DG>& dg, ModelArray& ma, const ModelArray::MultiDim& sliceIndices)
    {
        // Use the 2D optimized version if possible
        if (ma.nDimensions() == 2) return dg2ma(dg, ma);

        return dg2ma(dg, ma, kFromLocation(sliceIndices, ma));
    }

    template <int N> static DGVector<N>& hField2dg(const HField& h, DGVector<N>& dg)
    {
        dg.col(0) = h.data();
        return dg;
    }

    template <int N> static HField& dg2hField(const DGVector<N>& dg, HField& h)
    {
        h.setData(dg.col(0));
        return h;
    }
private:
    static size_t kFromLocation(const ModelArray::MultiDim& loc, const ModelArray& ma)
    {
        const size_t twoD = 2;

        // Copy the slice index and prefix with two zeroes
        ModelArray::MultiDim fullDims(loc);
        for (size_t dim = 0; dim < twoD; ++dim) {
            fullDims.insert(fullDims.begin(), 0);
        }
        // Get the index of the first element of the slice
        size_t kIndex = ma.indexFromLocation(fullDims);
        // Divide by the size of the 2D slice (first two dimensions)
        for (size_t dim = 0; dim < twoD; ++dim) {
            kIndex /= ma.dimensions()[dim];
        }
        return kIndex;
    }

};

template <> inline DGVector<1>& DGModelArray::hField2dg(const HField& h, DGVector<1>& dg)
{
    return ma2dg(h, dg);
}

template <> inline HField& DGModelArray::dg2hField(const DGVector<1>& dg, HField& h)
{
    return dg2ma(dg, h);
}
} /* namespace Nextsim */

#endif /* DGMODELARRAY_HPP */
