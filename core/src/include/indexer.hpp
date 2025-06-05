/*!
 * @file indexer.hpp
 *
 * @date Nov 5, 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef INDEXER_HPP
#define INDEXER_HPP

#include <cstddef>
#include <limits>
#include <vector>

namespace Indexer {

/*!
 * @brief Produces a one-dimensional index of an array given a multi-dimensional location.
 *
 * @detail Given the dimensions of a mulit-dimensional array and a vector of positions across those
 * dimensions, this function returns the corresponding index of that location within the flattened,
 * one-dimensional representation of the array, as it would typically be stored in memory, or when
 * a multi-dimensional spatial array is flattened into a single logical dimension.
 *
 * @param dims An ordered container holding the sequence of dimension sizes.
 * @param loc An ordered container holding the position in each dimension.
 * @return The flattened one-dimensional index corresponding to the two input arguments.
 */
template <typename T = size_t, typename D = std::vector<T>, typename L = D>
T indexer(const D& dims, const L& loc)
{
    auto iloc = std::begin(loc);
    auto idim = std::begin(dims);
    T stride = 1;
    T index = 0;
    do {
        index += stride * *iloc;
        stride *= *idim;

    } while (++iloc != std::end(loc) && ++idim != std::end(dims));
    return index;
}

/*!
 * @brief Produces a location in a multi-dimensional array from a one-dimensional index.
 *
 * @details Given the dimensions of a multi-dimensional array and a one-dimensional index in the
 * flattened representation of that array, this function returns the equivalent location as a list
 * of positions on the axes of the array.
 *
 * @param dims An ordered container holding the sequence of dimension sizes.
 * @param index The one-dimensional index to be converted to a mulit-dimensional location.
 * @return The multi-dimensional location equivalent to the one-dimensional index with an array of
 *         the provided arguments.
 */
template <typename T = size_t, typename D = std::vector<T>, typename L = D>
L deIndexer(const D& dims, const T& index)
{
    L loc;
    loc.resize(dims.size());
    auto locIter = loc.begin();
    T runningIndex = index;
    for (auto dim : dims) {
        *locIter = runningIndex % dim;
        runningIndex /= dim;
        ++locIter;
    }
    return loc;
}
} /* namespace Indexer */
#endif /* INDEXER_HPP */
