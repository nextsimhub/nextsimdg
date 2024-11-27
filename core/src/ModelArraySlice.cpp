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

    copySliceWithIters(other.data, otherIter, data, thisIter);
    return *this;
}

ModelArraySlice& ModelArraySlice::operator=(ModelArray& ma)
{
    copyBetweenMAandMASlice(ma, *this, true, "ModelArraySlice::operator=(ModelArray)");
    return *this;
}

ModelArray& ModelArraySlice::copyToModelArray(ModelArray& target) const
{
    copyBetweenMAandMASlice(target, *this, false, "ModelArraySlice::copyToModelArray(ModelArray)");
    return target;
}

ModelArray::DataType& ModelArraySlice::copyToDataSlice(
    ModelArray::DataType& target, SliceIter& targetIter) const
{
    SliceIter iter(slice, data.dimensions());

    copySliceWithItersData(data.m_data, iter, target, targetIter);
    return target;
}

ModelArraySlice& ModelArraySlice::copyFromDataSlice(
    const ModelArray::DataType& source, SliceIter& sourceIter)
{
    SliceIter iter(slice, data.dimensions());

    copySliceWithItersData(source, sourceIter, data.m_data, iter);
    return *this;
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
        copySliceWithIters(ma, maIter, mas.data, masIter);
    } else {
        copySliceWithIters(mas.data, masIter, ma, maIter);
    }
}

void ModelArraySlice::copySliceWithIters(
    ModelArray& source, SliceIter& sourceIter, ModelArray& target, SliceIter targetIter)
{
    copySliceWithItersData(source.m_data, sourceIter, target.m_data, targetIter);
}

void ModelArraySlice::copySliceWithItersData(const ModelArray::DataType& source,
    SliceIter& sourceIter, ModelArray::DataType& target, SliceIter targetIter)
{
    const size_t targetNEl = targetIter.nElements(0);
    const size_t sourceNEl = sourceIter.nElements(0);

    while (!targetIter.isEnd() && !sourceIter.isEnd()) {
        const size_t targetIndex = targetIter.index();
        const size_t sourceIndex = sourceIter.index();
        target(Eigen::seqN(targetIndex, targetNEl, targetIter.step(0)), Eigen::all)
            = source(Eigen::seqN(sourceIndex, sourceNEl, sourceIter.step(0)), Eigen::all);
        targetIter.incrementDim(1);
        sourceIter.incrementDim(1);
    }
}

MASIter::MASIter(ModelArraySlice& mas)
    : data(mas.data)
    , iter(mas.slice, mas.data.dimensions())
{
}
} // namespace Nextsim
