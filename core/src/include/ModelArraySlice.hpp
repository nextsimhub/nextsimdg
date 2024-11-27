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

// Inheriting from std::iterator is deprecated after C++17
class MASIter {
public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = double;
    using difference_type = std::ptrdiff_t;
    using pointer = double*;
    using reference = double&;
    using Slice = ArraySlicer::Slice;
    using SliceIter = ArraySlicer::SliceIter;

    MASIter(ModelArraySlice& mas);
    double& operator*() const { return data[iter.index()]; }
    MASIter& operator++()
    {
        ++iter;
        return *this;
    }
    bool operator!=(const MASIter& other) const
    {
        //        std::cout << this->iter << "!=" << other.iter << "?" << std::endl;
        // Account for different end iters possibly not equating
        bool iterEquality = (iter == other.iter) || (iter.isEnd() && other.iter.isEnd());
        return (&data != &other.data) || !iterEquality;
    }
    bool operator==(const MASIter& other) const { return !(*this != other); }
    double* operator->() const { return &data[iter.index()]; }
    MASIter operator++(int)
    {
        MASIter copy = *this;
        ++(*this);
        return copy;
    }

    std::ostream& print(std::ostream& os) const { return os << "{" << &data << ":" << iter << "}"; }
    friend ModelArraySlice;

private:
    ModelArray& data;
    SliceIter iter;
};

inline std::ostream& operator<<(std::ostream& os, const MASIter& it) { return it.print(os); }
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

    // Assign data to this slice from a ModelArray
    ModelArraySlice& operator=(ModelArray& ma);
    // Copy data from another ModelArraySlice
    ModelArraySlice& operator=(ModelArraySlice& other);
    // Assign a scalar to the entire slice
    ModelArraySlice& operator=(double v);
    // Assign the contents of a buffer to a slice
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

    ModelArray& copyToModelArray(ModelArray& target) const;
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

    ModelArray::DataType& copyToDataSlice(
        ModelArray::DataType& target, SliceIter& targetIter) const;
    ModelArraySlice& copyFromDataSlice(const ModelArray::DataType& source, SliceIter& sourceIter);

    iterator begin();
    iterator end();

    Slice slice;

    friend MASIter;

private:
    static void copyBetweenMAandMASlice(
        ModelArray& ma, const ModelArraySlice& mas, bool toSlice, const std::string& functionName);
    static void copySliceWithIters(
        ModelArray& source, SliceIter& sourceIter, ModelArray& target, SliceIter targetIter);
    static void copySliceWithItersData(const ModelArray::DataType& source, SliceIter& sourceIter,
        ModelArray::DataType& target, SliceIter targetIter);

    ModelArray& data;
};

} // namespace Nextsim

#endif /* MODELARRAYSLICE_HPP */
