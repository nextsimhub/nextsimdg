/*!
 * @file EigenSlice.hpp
 *
 * @date 16 Jul 2025
 * @author Tom Meltzer <tdm39@cam.ac.uk>
 */

#ifndef EIGENSLICE_HPP
#define EIGENSLICE_HPP

#include "include/ModelArray.hpp"
#include "include/ModelArraySlice.hpp"
#include "include/Slice.hpp"
#include "include/dgVector.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>

#ifndef DGCOMP
#define DGCOMP 3
#endif

using Slice = ArraySlicer::Slice;
using SliceIter = ArraySlicer::SliceIter;

namespace Nextsim {

/**
 * @class EigenSlice
 * @brief Provides a view and slicing operations for Eigen-based data structures.
 *
 * This class enables slicing and manipulation of Eigen matrix-like or array-like objects,
 * such as those used in DGVector and ModelArray. It supports copying between slices,
 * printing reshaped data, and provides an interface for working with multi-dimensional
 * data in a flexible way.
 *
 * @tparam Derived The Eigen type derived from Eigen::DenseBase.
 *
 * @note The class is designed to work with both matrix-like and array-like Eigen objects,
 * depending on the instantiation.
 */
template <typename Derived> class EigenSlice {
private:
    // Depending on instatiation this can either be
    // * matrix-like object for DGVector
    // * array-like object for ModelArray
    Eigen::DenseBase<Derived>& data;

    Slice slice;
    using MultiDim = std::vector<size_t>;

    // spatial dimensions and DoF
    const MultiDim dimensions;

public:
    explicit EigenSlice(ModelArray& ma, const Slice& sl)
        : data(ma.m_data)
        , dimensions(ma.dimensions())
        , slice(sl)
    {
    }

    explicit EigenSlice(DGVector<DGCOMP>& dgv, const Slice& sl, const ParametricMesh& smesh)
        : data(dgv)
        , dimensions({ smesh.nx, smesh.ny })
        , slice(sl)
    {
    }

    /*!
     * @brief Copies data from EigenSlice to target buffer
     *
     * @param target The target buffer
     * @return A reference to the original EigenSlice object.
     */
    template <typename S> const EigenSlice& copyToSlicedBuffer(EigenSlice<S>& target) const
    {
        SliceIter iter(slice, dimensions);
        SliceIter targetIter(target.slice, target.dimensions);

        ModelArraySlice::copySliceWithIters(data, iter, target.data, targetIter);
        // Return this, even though it is unchanged. The code looks weird if
        // the return value is the target buffer.
        return *this;
    }

    /*!
     * @brief Copies data to EigenSlice from source buffer
     *
     * @param source The source buffer
     * @return A reference to the updated EigenSlice object.
     */
    template <typename S> EigenSlice& copyFromSlicedBuffer(const EigenSlice<S>& source)
    {
        SliceIter iter(slice, dimensions);
        SliceIter sourceIter(source.slice, source.dimensions);

        ModelArraySlice::copySliceWithIters(source.data, sourceIter, data, iter);
        return *this;
    }

    /**
     * @brief Prints each column of the data as a reshaped 2D matrix.
     *
     * Iterates over all columns of the internal Eigen data object, reshaping each column
     * into a matrix of size (dimX, dimY) as specified by the dimensions member.
     * Each reshaped column is printed to the standard output, prefixed by its DG element index.
     */
    void print()
    {
        // Assuming data is a 2D Eigen object and dimensions = { rows, cols }
        size_t dimX = dimensions[0];
        size_t dimY = dimensions[1];

        // Reshape and print each column
        for (size_t col = 0; col < data.cols(); col++) {
            std::cout << "DG element " << col << std::endl;
            std::cout << data.col(col).reshaped(dimX, dimY) << std::endl;
        }
    }

private:
    // We need to add different Templated versions explicitly as friend
    // to allow access to private variables, as intended.
    template <typename D> friend class EigenSlice;
};

} // namespace Nextsim

#endif /* EIGENSLICE_HPP */
