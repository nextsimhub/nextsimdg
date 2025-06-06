/*!
 * @file ModelArraySlice.cpp
 *
 * @date Nov 11, 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ModelArraySlice.hpp"

namespace Nextsim {

ModelArraySlice& ModelArraySlice::operator=(double v)
{
    /*
     * The left most index varies fastest, as is the ModelArray standard
     * outside of components. This can be directly used as arguments to an
     * Eigen::seq object
     */
    SliceIter si(slice, data.dimensions());
    while (!si.isEnd()) {
        const size_t index = si.index();
        data.m_data(Eigen::seqN(index, si.nElements(0), si.step(0)), Eigen::all) = v;
        si.incrementDim(1);
    }
    return *this;
}

ModelArraySlice& ModelArraySlice::operator=(ModelArraySlice& other)
{
    SliceIter thisIter(slice, data.dimensions());
    SliceIter otherIter(other.slice, other.data.dimensions());
    // Check that the shapes match
    auto thisShape = thisIter.shape();
    auto otherShape = otherIter.shape();
    if (thisShape.size() != otherShape.size())
        throw std::out_of_range("ModelArraySlice dimension mismatch");
    for (auto i = 0; i < thisShape.size(); ++i) {
        if (thisShape[i] != otherShape[i])
            throw std::out_of_range("ModelArraySlice shape mismatch");
    }

    copySliceWithIters(other.data.m_data, otherIter, data.m_data, thisIter);
    return *this;
}

ModelArraySlice& ModelArraySlice::operator=(ModelArray& ma)
{
    copyBetweenMAandMASlice(ma, *this, true, "ModelArraySlice::operator=(ModelArray)");
    return *this;
}

ModelArraySlice& ModelArraySlice::operator=(const ModelArray::DataType& buffer)
{
    SliceIter thisIter(slice, data.dimensions());
    // The iterator for the buffer must have the same dimensions as the target slice.
    Slice::VBounds bufferBounds(data.nDimensions(), Slice::Bounds());
    ModelArray::MultiDim bufferDims(data.nDimensions());
    for (size_t dim = 0; dim < data.nDimensions(); ++dim) {
        bufferDims[dim] = thisIter.nElements(dim);
    }
    SliceIter bufferIter(bufferBounds, bufferDims);
    copySliceWithIters(buffer, bufferIter, data.m_data, thisIter);

    return *this;
}

ModelArraySlice::operator ModelArray::DataType()
{
    SliceIter thisIter(slice, data.dimensions());
    Slice::VBounds bufferBounds(data.nDimensions(), Slice::Bounds());
    ModelArray::MultiDim bufferDims(data.nDimensions());
    size_t nElements = 1;
    for (size_t dim = 0; dim < data.nDimensions(); ++dim) {
        bufferDims[dim] = thisIter.nElements(dim);
        nElements *= bufferDims[dim];
    }
    SliceIter bufferIter(bufferBounds, bufferDims);
    ModelArray::DataType buffer;
    buffer.resize(nElements, data.nComponents());
    copySliceWithIters(data.m_data, thisIter, buffer, bufferIter);

    return buffer;
}

ModelArray& ModelArraySlice::copyToModelArray(ModelArray& target) const
{
    copyBetweenMAandMASlice(target, *this, false, "ModelArraySlice::copyToModelArray(ModelArray)");
    return target;
}

ModelArraySlice::iterator ModelArraySlice::begin()
{
    iterator iter(*this);
    iter.iter.toBegin(); // To ensure we are at the beginning, in case SliceIter is ever modified.
    return iter;
}

ModelArraySlice::iterator ModelArraySlice::end()
{
    iterator iter(*this);
    iter.iter.toEnd();
    return iter;
}

void ModelArraySlice::copyBetweenMAandMASlice(
    ModelArray& ma, const ModelArraySlice& mas, bool toSlice, const std::string& functionName)
{
    // Shape check
    SliceIter masIter(mas.slice, mas.data.dimensions());

    // The slice needs to have at least as many dimensions as the ModelArray
    for (size_t dim = 0; dim < ma.nDimensions(); ++dim) {
        if (masIter.nElements(dim) != ma.dimensions()[dim])
            throw std::domain_error(functionName + ": shape mismatch.");
    }
    // But any extra dimensions must have length 1
    if (mas.slice.n() > ma.nDimensions()) {
        for (size_t dim = ma.nDimensions(); dim < mas.slice.n(); ++dim) {
            if (masIter.nElements(dim) != 1)
                throw std::domain_error(
                    functionName + ": additional dimensions must have length 1.");
        }
    }

    Slice::Bounds all = { {}, {} };
    Slice::Bounds first = { 0 };
    Slice::VBounds bounds;
    bounds.resize(mas.slice.n());
    SliceIter::MultiDim extendedDims;
    extendedDims.resize(mas.slice.n());
    for (size_t dim = 0; dim < ma.nDimensions(); ++dim) {
        bounds[dim] = all;
        extendedDims[dim] = ma.dimensions()[dim];
    }
    for (size_t dim = ma.nDimensions(); dim < mas.slice.n(); ++dim) {
        bounds[dim] = first;
        extendedDims[dim] = 1;
    }

    Slice maSlice(bounds);
    SliceIter maIter(maSlice, extendedDims);

    // Copy the data in the correct direction
    if (toSlice) {
        copySliceWithIters(ma.m_data, maIter, mas.data.m_data, masIter);
    } else {
        copySliceWithIters(mas.data.m_data, masIter, ma.m_data, maIter);
    }
}

MASIter::MASIter(ModelArraySlice& mas)
    : data(mas.data)
    , iter(mas.slice, mas.data.dimensions())
{
}
} // namespace Nextsim
