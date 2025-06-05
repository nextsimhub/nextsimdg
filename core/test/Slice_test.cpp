/*!
 * @file Slice_test.cpp
 *
 * @date 5 Nov 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "include/Slice.hpp"

#include <array>

#define ceil(num , denom) (((num) + (denom) - 1) / (denom))

using Slice = ArraySlicer::Slice;
using SliceIter = ArraySlicer::SliceIter;
using MultiDim = SliceIter::MultiDim;

std::string to_string(const MultiDim& v)
{
    std::string s = "(";
    for (auto dim = 0; dim < v.size(); ++dim) {
        s += std::to_string(v[dim]);
        if (dim != v.size() - 1) s += ",";
    }
    return s += ")";
}

std::string to_string(const Slice& l)
{
    std::string s = "(";
    // Use the "[" is inclusive, ")" is exclusive notation.
    for (auto dim = 0; dim < l.n(); ++dim) {
        s += "[" + (l.bounds[dim].start.isAll() ? "" : std::to_string(l.bounds[dim].start))
                + ":"+ (l.bounds[dim].stop.isAll() ? "" : std::to_string(l.bounds[dim].stop)) + ")";
        if (dim != l.n() - 1) s += ",";
    }
    return s += ")";
}

TEST_SUITE_BEGIN("Indexer");
TEST_CASE("indexer <-> deIndexer")
{
    using MultiDim = std::vector<size_t>;
    MultiDim origLoc = { 2, 3, 5, 7, 11, 13 };
    MultiDim dims = { 4, 6, 8, 10, 12, 14 };

    size_t index = Indexer::indexer(dims, origLoc);
    MultiDim finalLoc = Indexer::deIndexer(dims, index);
    for (size_t i = 0; i < dims.size(); ++i) {
        REQUIRE(origLoc[i] == finalLoc[i]);
    }
}
TEST_SUITE_END();

TEST_SUITE_BEGIN("Slice");
TEST_CASE("One dimensional indexing")
{
    // A single element
    Slice element3 {{{3U}}};
    SliceIter iter3(element3, {10});
    REQUIRE(iter3.index() == 3);
    REQUIRE(iter3.size() == 1);

    // All the elements of an array
    // Implicitly
    SliceIter impAll({{{}}}, {3});
    size_t count = 0;
    while(!impAll.isEnd()) {
        ++impAll;
        ++count;
    }
    REQUIRE(count == 3);
    // Explicitly
    SliceIter expAll({{{{}}}}, {3});
    count = 0;
    while(!expAll.isEnd()) {
        ++expAll;
        ++count;
    }
    REQUIRE(count == 3);

    // Step backwards
    SliceIter expAllB({{{{}, {}, -1}}}, {3});
    count = 0;
    while(!expAllB.isEnd()) {
        ++expAllB;
        ++count;
    }
    REQUIRE(count == 3);

    // A range
    Slice range36 {{{3, 6}}};
    SliceIter iter36(range36, {10});
    REQUIRE(iter36.index() == 3);
    REQUIRE(iter36.isBegin());
    ++iter36;
    REQUIRE(iter36.index() == 4);
    // post incrementing
    REQUIRE(iter36++.index() == 4);
    REQUIRE(iter36.index() == 5);
    ++iter36;
    REQUIRE(iter36.isEnd());
    // Increment past the end
    ++iter36;

    // range to the end of the array
    Slice range3_ {{{3, {}}}};
    SliceIter iter3_(range3_, {10});
    // count the values until the end 3456789
    count = 0;
    while(!iter3_.isEnd()) {
        ++iter3_;
        ++count;
        REQUIRE(count < 10);
    }
    REQUIRE(iter3_.isEnd());
    REQUIRE(iter3_.size() == 7);
    REQUIRE(count == 7);

    // Non-unit stride, also full length of the array.
    Slice stride2 {{{{},{},2}}};
    SliceIter iter__2(stride2, {10});
    REQUIRE(iter__2.index() == 0);
    count = 0;
    size_t expt = 10 / 2;
    REQUIRE(iter__2.size() == expt);
    while (!iter__2.isEnd()) {
        REQUIRE(iter__2.index() % 2 == 0);
        ++iter__2;
        ++count;
    }
    REQUIRE(count == expt);

    // And starting from 1
    Slice stride2a = {{{1, {}, 2}}};
    SliceIter iter1_2(stride2a, {10});
    REQUIRE(iter1_2.index() == 1);
    count = 0;
    // expt still = 5
    while (!iter1_2.isEnd()) {
        REQUIRE(iter1_2.index() % 2 == 1);
        ++iter1_2;
        ++count;
    }
    REQUIRE(count == expt);

    Slice allOneD;
}

TEST_CASE("Multidimensional indexing")
{
    // A single element
    Slice element {{{3}, {5}}};
    SliceIter iter610(element, {6, 10});
    REQUIRE(iter610.index() == Indexer::indexer({6, 10}, {3, 5}));
    REQUIRE(SliceIter(element, {8, 7}).index() == Indexer::indexer({8, 7}, {3, 5}));

    // Check that a mismatch in number of dimensions is correctly detected
    REQUIRE_THROWS_AS(SliceIter(element, {8, 7, 6}), std::invalid_argument);

    // A multidimensional slice
    Slice elements3748 {{{3, 7}, {4, 8}}};
    const size_t nx = 11;
    const size_t ny = 13;
    SliceIter iter3748(elements3748, {nx, ny});
    size_t index = iter3748.index();
    REQUIRE(index == Indexer::indexer({nx, ny}, {3, 4}));
    size_t count = 0;
    size_t indexLast = index - 1;
    while(!iter3748.isEnd()) {
        index = iter3748.index();
        std::ptrdiff_t deltaIndex = index - indexLast;
        bool xStep = deltaIndex == 1;
        bool yStep = deltaIndex == (Indexer::indexer({nx, ny}, {3, 5}) - Indexer::indexer({nx, ny}, {6, 4}));
        bool allowedStep = xStep || yStep;
        REQUIRE_MESSAGE(allowedStep, "Forbidden step value Δi=", deltaIndex);

        REQUIRE(iter3748.position().size() == 2);
        REQUIRE(iter3748.position()[0] >= 3);
        REQUIRE(iter3748.position()[0] < 7);
        REQUIRE(iter3748.position()[1] >= 4);
        REQUIRE(iter3748.position()[1] < 8);

        REQUIRE(count <= 16);
        indexLast = index;
        ++iter3748;
        ++count;
    }
    REQUIRE(iter3748.size() == 16);

    // 8 dimensional array slicing
    Slice elements8d {{{2, 6}, {4, 9}, {6, 12}, {8, 15}, {10, 18}, {12, 21}, {14, 24}, {16, 27}}};
    std::vector<size_t> ni = { 7, 12, 13, 16, 30, 30, 30, 30};
    SliceIter iter8d(elements8d, ni);
    count = 0;
    const size_t expt = 11 * 10 * 9 * 8 * 7 * 6 * 5 * 4;
    REQUIRE(iter8d.size() == expt);
    std::string message = "";
    while (!iter8d.isEnd()) {
        auto loc = Indexer::deIndexer(ni, iter8d.index());
        for (size_t dim = 0; dim < ni.size(); ++dim) {
             if (!(loc[dim] >= elements8d.bounds[dim].start)) {
                 message = "start";
             } else if (!(loc[dim] < elements8d.bounds[dim].stop)) {
                 message = "stop";
             }
             if (!message.empty()) goto end8d;
        }
        ++iter8d;
        ++count;
    }
end8d:
    REQUIRE_MESSAGE(message.empty(), "Position ", to_string(Indexer::deIndexer(ni, iter8d.index())), " is out of bounds ", to_string(elements8d));
    REQUIRE(count == expt);
    SliceIter::MultiDim targ = {4, 5, 6, 7, 8, 9, 10, 11};
    REQUIRE(iter8d.shape().size() == targ.size());
    for (auto i = 0; i < targ.size(); ++i) {
        REQUIRE(iter8d.shape()[i] == targ[i]);
    }

    // n to whatever in 2 dimensions
    const size_t xi = 3;
    const size_t yi = 5;
    Slice slice3_5_ {{{xi, {}}, {yi, {}}}};
    std::vector<std::vector<size_t>> dims = { { 11, 17 }, { 37, 43 } };
    for (std::vector<size_t> nj : dims) {
        count = 0;
        for (SliceIter iter(slice3_5_, nj); !iter.isEnd(); ++iter) {
            ++count;
        }
        REQUIRE(count == (nj[0] - xi) * (nj[1] - yi));
    }

    // Multiple dimensions, multiple strides
    auto dx = 2;
    auto dy = 3;
    auto lenY = 13;
    Slice sliceMultiStride {{{xi, {}, dx}, {yi, yi + lenY, dy}}};
    count = 0;
    auto dim = dims[1]; // Take the arrays sizes from above
    for (SliceIter iter(sliceMultiStride, dim); !iter.isEnd(); ++iter) {
        count++;
    }
    SliceIter::MultiDim shape = SliceIter(sliceMultiStride, dim).shape();
    SliceIter::MultiDim shapeTarget = {static_cast<size_t>(ceil(dim[0] - xi, dx)), static_cast<size_t>(ceil(lenY, dy))};
    REQUIRE(count == shapeTarget[0] * shapeTarget[1]);
    REQUIRE(shape.size() == shapeTarget.size());
    for (auto i = 0; i < shape.size(); ++i) {
        REQUIRE(shape[i] == shapeTarget[i]);
    }

    // Reuse elements8d to test higher dimensional incrementing
    size_t i1, i2;
    iter8d.toBegin();
    i1 = iter8d.index();
    i2 = (++iter8d).index();
    REQUIRE(Indexer::deIndexer(ni, i2)[0] == Indexer::deIndexer(ni, i1)[0] + 1);
    i1 = iter8d.incrementDim(1).index();
    REQUIRE(Indexer::deIndexer(ni, i1)[1] == Indexer::deIndexer(ni, i2)[1] + 1);
    i2 = iter8d.incrementDim(7).index();
    REQUIRE(Indexer::deIndexer(ni, i2)[7] == Indexer::deIndexer(ni, i1)[7] + 1);

    // Extract a single row of a 2d array
    const size_t nx1 = 12;
    const size_t ny1 = 15;
    SliceIter::MultiDim dims1 = {12, 15};
    Slice singleColumn {{{}, {0}}};
    count = 0;
    for (SliceIter iter(singleColumn, {nx1, ny1}); !iter.isEnd(); ++iter) {
        REQUIRE(Indexer::deIndexer(dims1, iter.index())[0] == count);
        ++count;
        REQUIRE(Indexer::deIndexer(dims1, iter.index())[1] == 0);
    }

}

TEST_CASE("Negative steps")
{
    Slice sdrawkcab {{{{}, {}, -1}}};
    size_t n = 8;
    SliceIter::MultiDim dim1d = {n};
    SliceIter reti(sdrawkcab, dim1d);
    // Use the start and nElements functions to get the real indices
    REQUIRE(reti.start() == n - 1);
    REQUIRE(reti.nElements() == n);
    REQUIRE(reti.step() == -1);

    Slice sdraw2 {{{6, 1, -2}}};
    SliceIter reti2(sdraw2, dim1d);
    REQUIRE(reti2.start() == 6);
    REQUIRE(reti2.nElements() == 3);
    REQUIRE(reti2.step() == -2);
    size_t count = 0;
    while(!reti2.isEnd()) {
        ++reti2;
        ++count;
    }
    REQUIRE(count == 3);
    REQUIRE(reti2.size() == 3);

    // Forwards in one dimension, backwards in another
    Slice to_fro {{{0, 8}, {7, {}, -1}}};
    MultiDim dims = {10, 10};
    SliceIter itereti(to_fro, dims);
    REQUIRE(itereti.index() == Indexer::indexer(dims, {0, 7}));
    REQUIRE(itereti.shape()[0] == 8);
    REQUIRE(itereti.shape()[1] == 8);
    count = 0;
    while(!itereti.isEnd()) {
        REQUIRE(Indexer::deIndexer(dims, itereti.index())[0] >= 0);
        REQUIRE(Indexer::deIndexer(dims, itereti.index())[0] < 8);
        REQUIRE(Indexer::deIndexer(dims, itereti.index())[1] >= 0);
        REQUIRE(Indexer::deIndexer(dims, itereti.index())[1] < 8);
        ++itereti;
        ++count;
    }
    REQUIRE(count == itereti.shape()[0] * itereti.shape()[1]);
    REQUIRE(itereti.size() == count);
}

bool matchIndices(SliceIter& si, std::vector<size_t> compare, bool debug = false)
{
    si.toBegin();
    size_t compareIndex = 0;
    while(!si.isEnd())
    {
        if (debug) std::cout << si.index() << "?="<< compare[compareIndex] << " ";
        if (compareIndex == compare.size()) return false;
        if (debug) std::cout << " √ compare.size";
        if (si.index() != compare[compareIndex]) return false;
        ++si;
        ++compareIndex;
        if (debug) std::cout << ". ";
    }
    if (debug) std::cout << std::endl;
    // Only return true if there are no more elements of compare remaining
    return (compareIndex == compare.size());
}

TEST_CASE("Test matchIndices")
{
    SliceIter allTen({{{{}, {}}}}, {10});
    REQUIRE(matchIndices(allTen, {0,1,2,3,4,5,6,7,8,9}));
    // Mismatched index
    REQUIRE_FALSE(matchIndices(allTen, {0,1,2,3,4,5,6000,7,8,9}));
    // Too few elements in compare
    REQUIRE_FALSE(matchIndices(allTen, {0,1,2,3,4,5,6,7,8}));
    // Too many elements in compare
    REQUIRE_FALSE(matchIndices(allTen, {0,1,2,3,4,5,6,7,8,9,9}));

    // Test backwards and a sub set
    SliceIter back({{{8, 3, -1}}}, {10});
    REQUIRE(matchIndices(back, {8, 7, 6, 5, 4}));
    // Multidimensional
    SliceIter twod({{{},{}}}, {2,2});
    SliceIter::MultiDim shape2 = twod.shape();
    REQUIRE(matchIndices(twod, {0,1,2,3}));
    SliceIter twodSlice({{{2, 4}, {3, 5}}}, {5, 6});
    REQUIRE(matchIndices(twodSlice, {17,18,22,23}));

    // Test an empty slice and an empty comparison vector
    SliceIter none({{{6, 3}}}, {10});
    REQUIRE(matchIndices(none, {}));
}

// Wrap the comparison boilerplate code
bool match12(const Slice::VBounds& b, std::vector<size_t> compare, bool debug = false)
{
    SliceIter a(b, {12});
    return matchIndices(a, compare, debug);
}

TEST_CASE("Indexing behaviour")
{
    // Compares the behaviour of the indexing to that of the numpy Python library.

    // step = 1
    // -ve start
    // stop < -dim
    REQUIRE(match12({{{-6, -16}}}, {}));
    // stop = -dim
    REQUIRE(match12({{{-6, -12}}}, {}));
    // -dim < stop < start
    REQUIRE(match12({{{-6, -8}}}, {}));
    // start < stop < 0
    REQUIRE(match12({{{-6, -4}}}, {6, 7}));
    // stop = 0
    REQUIRE(match12({{{-6, 0}}}, {}));
    // stop > 0, before start
    REQUIRE(match12({{{-6, 4}}}, {}));
    // stop > 0, after start
    REQUIRE(match12({{{-6, 8}}}, {6, 7}));
    // stop = length
    REQUIRE(match12({{{-6, 12}}}, {6, 7, 8, 9, 10, 11}));
    // stop < length
    REQUIRE(match12({{{-6, 16}}}, {6, 7, 8, 9, 10, 11}));
    // default stop
    REQUIRE(match12({{{-6, {}}}}, {6, 7, 8, 9, 10, 11}));
    // stop not set
    REQUIRE(match12({{{-6}}}, {6}));

    // 0 start
    // stop < -dim
    REQUIRE(match12({{{0, -16}}}, {}));
    // stop = -dim
    REQUIRE(match12({{{0, -12}}}, {}));
    // -dim < stop < start
    REQUIRE(match12({{{0, -8}}}, {0, 1, 2, 3}));
    // stop = 0
    REQUIRE(match12({{{0, 0}}}, {}));
    // 0 < stop, after start
    REQUIRE(match12({{{0, 8}}}, {0, 1, 2, 3, 4, 5, 6, 7}));
    // stop = length
    REQUIRE(match12({{{0, 12}}}, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}));
    // stop < length
    REQUIRE(match12({{{0, 16}}}, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}));
    // default stop
    REQUIRE(match12({{{0, {}}}}, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}));
    // stop not set
    REQUIRE(match12({{{0}}}, {0}));

    // +ve start
    // stop < -dim
    REQUIRE(match12({{{6, -16}}}, {}));
    // stop = -dim
    REQUIRE(match12({{{6, -12}}}, {}));
    // stop < 0, before start
    REQUIRE(match12({{{6, -8}}}, {}));
    // stop < 0, after start
    REQUIRE(match12({{{6, -4}}}, {6, 7}));
    // stop = 0
    REQUIRE(match12({{{6, 0}}}, {}));
    // 0 < stop < start
    REQUIRE(match12({{{6, 4}}}, {}));
    // stop > start
    REQUIRE(match12({{{6, 8}}}, {6, 7}));
    // stop = length
    REQUIRE(match12({{{6, 12}}}, {6, 7, 8, 9, 10, 11}));
    // stop < length
    REQUIRE(match12({{{6, 16}}}, {6, 7, 8, 9, 10, 11}));
    // default stop
    REQUIRE(match12({{{6, {}}}}, {6, 7, 8, 9, 10, 11}));
    // stop not set
    REQUIRE(match12({{{6}}}, {6}));

    // default start
    // stop < -dim
    REQUIRE(match12({{{{}, -16}}}, {}));
    // stop = -dim
    REQUIRE(match12({{{{}, -12}}}, {}));
    // stop < 0
    REQUIRE(match12({{{{}, -8}}}, {0, 1, 2, 3}));
    // stop = 0
    REQUIRE(match12({{{{}, 0}}}, {}));
    // stop > 0
    REQUIRE(match12({{{{}, 8}}}, {0, 1, 2, 3, 4, 5, 6, 7}));
    // stop = length
    REQUIRE(match12({{{{}, 12}}}, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}));
    // stop < length
    REQUIRE(match12({{{{}, 16}}}, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}));
    // default stop
    REQUIRE(match12({{{{}, {}}}}, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}));
    // stop not set
    SliceIter def({{{{}}}}, {12});
    REQUIRE(matchIndices(def, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}));
    // TODO get this test working with match12
//    REQUIRE(match12(Slice::VBounds({{{{}}}}), {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}, true));

    // step = -1
    // -ve start
    // stop < -dim
    REQUIRE(match12({{{-6, -16, -1}}}, {6, 5, 4, 3, 2, 1, 0}));
    // stop = -dim
    REQUIRE(match12({{{-6, -12, -1}}}, {6, 5, 4, 3, 2, 1}));
    // -dim < stop < start
    REQUIRE(match12({{{-6, -8, -1}}}, {6, 5}));
    // start < stop < 0
    REQUIRE(match12({{{-6, -4, -1}}}, {}));
    // stop = 0
    REQUIRE(match12({{{-6, 0, -1}}}, {6, 5, 4, 3, 2, 1}));
    // stop > 0, before start
    REQUIRE(match12({{{-6, 4, -1}}}, {6, 5}));
    // stop > 0, after start
    REQUIRE(match12({{{-6, 8, -1}}}, {}));
    // stop = length
    REQUIRE(match12({{{-6, 12, -1}}}, {}));
    // stop < length
    REQUIRE(match12({{{-6, 16, -1}}}, {}));
    // default stop
    REQUIRE(match12({{{-6, {}, -1}}}, {6, 5, 4, 3, 2, 1, 0}));

    // 0 start
    // stop < -dim
    REQUIRE(match12({{{0, -16, -1}}}, {0}));
    // stop = -dim
    REQUIRE(match12({{{0, -12, -1}}}, {}));
    // -dim < stop < start
    REQUIRE(match12({{{0, -8, -1}}}, {}));
    // stop = 0
    REQUIRE(match12({{{0, 0, -1}}}, {}));
    // 0 < stop, after start
    REQUIRE(match12({{{0, 8, -1}}}, {}));
    // stop = length
    REQUIRE(match12({{{0, 12, -1}}}, {}));
    // stop < length
    REQUIRE(match12({{{0, 16, -1}}}, {}));
    // default stop
    REQUIRE(match12({{{0, {}, -1}}}, {0}));

    // +ve start
    // stop < -dim
    REQUIRE(match12({{{6, -16, -1}}}, {6, 5, 4, 3, 2, 1, 0}));
    // stop = -dim
    REQUIRE(match12({{{6, -12, -1}}}, {6, 5, 4, 3, 2, 1}));
    // stop < 0, before start
    REQUIRE(match12({{{6, -8, -1}}}, {6, 5}));
    // stop < 0, after start
    REQUIRE(match12({{{6, -4, -1}}}, {}));
    // stop = 0
    REQUIRE(match12({{{6, 0, -1}}}, {6, 5, 4, 3, 2, 1}));
    // 0 < stop < start
    REQUIRE(match12({{{6, 4, -1}}}, {6, 5}));
    // stop > start
    REQUIRE(match12({{{6, 8, -1}}}, {}));
    // stop = length
    REQUIRE(match12({{{6, 12, -1}}}, {}));
    // stop < length
    REQUIRE(match12({{{6, 16, -1}}}, {}));
    // default stop
    REQUIRE(match12({{{6, {}, -1}}}, {6, 5, 4, 3, 2, 1, 0}));

    // default start
    // stop < -dim
    REQUIRE(match12({{{{}, -16, -1}}}, {11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0}));
    // stop = -dim
    REQUIRE(match12({{{{}, -12, -1}}}, {11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1}));
    // stop < 0
    REQUIRE(match12({{{{}, -8, -1}}}, {11, 10, 9, 8, 7, 6, 5}));
    // stop = 0
    REQUIRE(match12({{{{}, 0, -1}}}, {11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1}));
    // stop > 0
    REQUIRE(match12({{{{}, 8, -1}}}, {11, 10, 9}));
    // stop = length
    REQUIRE(match12({{{{}, 12, -1}}}, {}));
    // stop < length
    REQUIRE(match12({{{{}, 16, -1}}}, {}));
    // default stop
    REQUIRE(match12({{{{}, {}, -1}}}, {11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0}));

}

// Creates a Slice with step = 0
bool makeZeroStep()
{
    Slice slice{{{{},{},0}}};
    return true;
}

TEST_CASE("Zero step")
{
    REQUIRE_THROWS(makeZeroStep());
}

TEST_CASE("SliceIter equality")
{
    MultiDim dims = {17, 19};
    SliceIter a({{{2, -2}, {3, -3}}}, dims);
    // Different number of dimensions
    REQUIRE_FALSE(a == SliceIter({{{2, -2}, {3, -3}, {}}}, {17, 19, 23}));
    // Different dimensions
    REQUIRE_FALSE(a == SliceIter({{{2, -2}, {3, -3}}}, {17, 18}));
    // Different limits
    REQUIRE_FALSE(a == SliceIter({{{2, -2}, {3, -2}}}, dims));
    // Different position
    SliceIter b({{{2, -2}, {3, -3}}}, dims);
    ++a;
    REQUIRE_FALSE(a == b);
    ++b;
    REQUIRE(a == b);

}

TEST_CASE("Two dimensions, slicing the end of the second dimension")
{
    MultiDim dims = {11, 13};
    Slice top {{{},{-1}}};
    SliceIter topIter(top, dims);
    topIter.toBegin();
    REQUIRE_FALSE(topIter.isEnd());
    REQUIRE(topIter.index() == Indexer::indexer(dims, {0, 12}));
    ++topIter;
    REQUIRE(topIter.index() == Indexer::indexer(dims, {1, 12}));
}
TEST_SUITE_END();
