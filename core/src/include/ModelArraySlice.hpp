/*!
 * @file ModelArraySlice.hpp
 *
 * @date Nov 8, 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef MODELARRAYSLICE_HPP
#define MODELARRAYSLICE_HPP

#include "include/ModelArray.hpp"
#include "include/Slice.hpp"

#include <iostream>
#include <iterator>

namespace Nextsim {
class MASIter;
class ModelArraySlice;
std::ostream& operator<<(std::ostream& os, const MASIter& it);

/*!
 * An iterator class for ModelArraySlices. Conforms to the forward iterator specifications.
 */
class MASIter {
    // Inheriting from std::iterator is deprecated after C++17
public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = double;
    using difference_type = std::ptrdiff_t;
    using pointer = double*;
    using reference = double&;
    using Slice = ArraySlicer::Slice;
    using SliceIter = ArraySlicer::SliceIter;

    MASIter(ModelArraySlice& mas);
    /*!
     * Provides read/write access to the data at the current location of the iterator.
     * @return reference access to the data the iterator refers to.
     */
    double& operator*() const { return data[iter.index()]; }
    //! Increments the position of the iterator.
    MASIter& operator++()
    {
        ++iter;
        return *this;
    }
    /*!
     * Checks if another iterator is no equal to this.
     * @param the second iterator to compare.
     * @return true when the two iterators are different.
     */
    bool operator!=(const MASIter& other) const
    {
        // Account for different end iters possibly not equating
        bool iterEquality = (iter == other.iter) || (iter.isEnd() && other.iter.isEnd());
        return (&data != &other.data) || !iterEquality;
    }
    /*!
     * Checks if another iterator is equal to this.
     * @param the second iterator to compare
     * @return true when the two iterators point to the same point of the same array.
     */
    bool operator==(const MASIter& other) const { return !(*this != other); }
    /*!
     * Provides member access for the data at the current location of the iterator.
     * @return a pointer to the current data
     */
    double* operator->() const { return &data[iter.index()]; }
    /*!
     * Post-increments the iterator.
     * @return The state of the iterator before it was incremented.
     */
    MASIter operator++(int)
    {
        MASIter copy = *this;
        ++(*this);
        return copy;
    }

    /*!
     * Pretty-prints the current state of the iterator, including the address of the data array.
     */
    std::ostream& print(std::ostream& os) const { return os << "{" << &data << ":" << iter << "}"; }
    friend ModelArraySlice;

private:
    ModelArray& data;
    SliceIter iter;
};

/*!
 * Inserts a printed representation of the iterator into a stream.
 */
inline std::ostream& operator<<(std::ostream& os, const MASIter& it) { return it.print(os); }
/*!
 * A class to provide slicing of ModelArray data arrays.
 */
class ModelArraySlice {
public:
    using Slice = ArraySlicer::Slice;
    using SliceIter = ArraySlicer::SliceIter;
    using iterator = MASIter;
    ModelArraySlice() = delete;
    ModelArraySlice(ModelArray& ma, const Slice& sl)
        : data(ma)
        , slice(sl)
    {
    }

    /*!
     * Assigns data to this slice from a ModelArray.
     * @param ma The ModelArray object to read the data from
     * @return A reference to the updated ModelArraySlice object
     */
    ModelArraySlice& operator=(ModelArray& ma);
    /*!
     * Copies data from another ModelArraySlice.
     * @param other The ModelArraySlice object to read the data from
     * @return A reference to the updated ModelArraySlice object.
     */
    ModelArraySlice& operator=(ModelArraySlice& other);
    /*!
     * Assigns a value to the entire slice, based on the given scalar value.
     * @param v The value to which to set the slice.
     * @return A reference to the updated ModelArraySlice.
     */
    ModelArraySlice& operator=(double v);
    /*!
     * Assigns the contents of a buffer to the slice. This function copies the entire buffer.
     * @param buffer The buffer to copy the data from. Can be of any type which provides a forward
     *               iterator. The provided object must have the same number of elements as the
     *               slice.
     * @return A reference to the updated ModelArraySlice object.
     */
    template <typename T> ModelArraySlice& operator=(const T& buffer)
    {
        // make no especial attempt at efficiency here
        SliceIter thisIter(slice, data.dimensions());
        auto biter = buffer.begin();

        while (!thisIter.isEnd()) {
            // If the buffer ends before the slice, throw an exception
            if (biter == buffer.end()) {
                throw std::length_error("ModelArraySlice::operator=(T): buffer exhausted");
            }
            data[thisIter.index()] = *biter;
            ++thisIter;
            ++biter;
        }
        return *this;
    }

    /*!
     * @brief Assigns to a slice from a instance of ModelArray::DataType.
     *
     * @details Fills the slice with data from an instance of ModelArray::DataType. The data source
     * should be sufficiently sized that the data is copied without being indexed out-of-bounds.
     * The number of components of the slice ModelArray and the size of the second dimension of the
     * source array should match. The spatial dimensions of the slice will copy from the first
     * dimension of the source array, treating it as a flattened one-dimensional version of the
     * slice, with the first spatial dimension of the slice varying fastest.
     *
     * @params dataBuffer The source of the data to be copied into the slice.
     */
    ModelArraySlice& operator=(const ModelArray::DataType& dataBuffer);

    /*!
     * @brief Returns the contents of this slice as an instance of ModelArray::DataType.
     *
     * @details Creates an instance of ModelArray::DataType with the contents of the slice,
     * including all components. Current ModelArray::DataType is Eigen::Array. In this case all
     * spatial dimensions are flattened to one dimension, as they are internally in ModelArray,
     * and occupy the first dimension of the Eigen::Array. Any components of the sliced array take
     * up the second dimension. That is, a 4 x 5 slice with 3 components will create a (20, 3)
     * Eigen::Array.
     */
    operator ModelArray::DataType();

    /*!
     * Copies the contents of a slice to a ModelArray with an equal number of elements.
     * @param target The ModelArray object to copy the contents of the slice to.
     * @return A reference to the updated ModelArray object.
     */
    ModelArray& copyToModelArray(ModelArray& target) const;
    /*!
     * Copies the contents of the slice to a buffer.
     * @param buffer The target of the copying. The buffer type must provide a forward operator.
     *               The buffer object must provide at least as many elements as exist in the
     *               slice.
     * @return A reference to the updated buffer object.
     */
    template <typename T> T& copyToBuffer(T& buffer)
    {
        // make no especial attempt at efficiency here
        SliceIter thisIter(slice, data.dimensions());
        auto biter = buffer.begin();

        while (!thisIter.isEnd()) {
            // If the buffer ends before the slice, throw an exception
            if (biter == buffer.end()) {
                throw std::length_error("ModelArraySlice::copyToBuffer(T): buffer exhausted");
            }
            *biter = data[thisIter.index()];
            ++thisIter;
            ++biter;
        }
        return buffer;
    }

private:
    // Copies data between any two objects which accept indexing by the Eigen seqN function. This
    // includes the array underlying ModelArray and the nextSIM dynamics array types.
    template <typename S, typename T = S>
    static void copySliceWithIters(
        const S& source, SliceIter& sourceIter, T& target, SliceIter targetIter)
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

public:
    /*!
     * Copies data from the slice to a slice of a buffer.
     * @param target The target buffer. The buffer type must be one that can be indexed by the
     *               Eigen::seqN function.
     * @param targetIter An instance of SliceIter created from the desired slice of the buffer and
     *                   the dimensions of the target object.
     * @return A reference to this ModelArraySlice (not the updated buffer).
     */
    template <typename T>
    const ModelArraySlice& copyToSlicedBuffer(T& target, SliceIter& targetIter) const
    {
        SliceIter iter(slice, data.dimensions());

        copySliceWithIters(data.m_data, iter, target, targetIter);
        // Return this, even though it is unchanged. The code looks weird if
        // the return value is the target buffer.
        return *this;
    }

    /*!
     * Copies data to the slice from a slice of a buffer.
     * @param source The source buffer.The buffer type must be one that can be indexed by the
     *               Eigen::seqN function.
     * @param sourceIter An instance of SliceIter created from the desired slice of the buffer and
     *                   the dimensions of the source object.
     * @return A reference to the updated ModelArraySlice object.
     */
    template <typename S>
    ModelArraySlice& copyFromSlicedBuffer(const S& source, SliceIter& sourceIter)
    {
        SliceIter iter(slice, data.dimensions());

        copySliceWithIters(source, sourceIter, data.m_data, iter);
        return *this;
    }

    /*!
     * Returns a forward iterator pointing to the start of the slice.
     */
    iterator begin();
    /*!
     * Returns a forward iterator pointing to one past the end of the slice.
     */
    iterator end();

    /*!
     * The slice this object represents. Updating this invalidates all current iterators.
     */
    Slice slice;

    friend MASIter;

private:
    static void copyBetweenMAandMASlice(
        ModelArray& ma, const ModelArraySlice& mas, bool toSlice, const std::string& functionName);

    ModelArray& data;
};

} // namespace Nextsim

#endif /* MODELARRAYSLICE_HPP */
