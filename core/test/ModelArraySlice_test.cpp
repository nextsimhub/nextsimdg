/*!
 * @file ModelArraySlice_test.cpp
 *
 * @date Nov 8, 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "include/ModelArraySlice.hpp"

#include <array>

using Slice = ArraySlicer::Slice;
using SliceIter = ArraySlicer::SliceIter;

const auto nx = 23;
const auto ny = 17;
const auto nz = 5;

namespace Nextsim {
TEST_SUITE_BEGIN("ModelArraySlice");
TEST_CASE("Assign scalar to slice")
{
    ModelArray::setDimension(ModelArray::Dimension::X, nx);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny);

    TwoDField target(ModelArray::Type::TWOD);
    target.resize();
    double orig = 0.;
    target = orig;
    Slice patch {{{5, 14}, {8, 15}}};
    ModelArraySlice mas(target, patch);
    double set = 1.;
    mas = set;
    REQUIRE(target(5, 8) == set);
    REQUIRE(target(13, 8) == set);
    REQUIRE(target(5, 14) == set);
    REQUIRE(target(13, 14) == set);

    REQUIRE(target(4, 8) == orig);
    REQUIRE(target(5, 7) == orig);
    REQUIRE(target(14, 8) == orig);
    REQUIRE(target(13, 7) == orig);
    REQUIRE(target(4, 14) == orig);
    REQUIRE(target(5, 15) == orig);
    REQUIRE(target(14, 14) == orig);
    REQUIRE(target(13, 15) == orig);
}

TEST_CASE("Slice to Slice")
{
    ModelArray::setDimension(ModelArray::Dimension::X, nx);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny);

    TwoDField source(ModelArray::Type::TWOD);
    source.resize();
    TwoDField target(ModelArray::Type::TWOD);
    target.resize();
    double targetv = -1.;
    target = targetv;

    for (size_t i = 0; i < source.size(); ++i) {
        source[i] = i;
    }

    // First test the error conditions
    ModelArraySlice masa(source, {{{1, 3}, {{}, 19}}});
    ModelArraySlice masb(source, {{{1, 3}, {8, 19}}});
    REQUIRE_THROWS_AS(masa = masb, std::out_of_range);

    size_t x0s = 0;
    size_t x1s = 3;
    size_t y0s = 0;
    size_t y1s = 14;

    size_t x0t = nx - x1s;
    size_t x1t = nx;
    size_t y0t = ny - y1s;
    size_t y1t = ny;
    ModelArraySlice sourceSlice(source, {{{x0s, x1s}, {y0s, y1s}}});
    ModelArraySlice targetSlice(target, {{{x0t, x1t}, {y0t, y1t}}});

    // Perform the Slice to Slice assignment
    targetSlice = sourceSlice;

    // Test the areas that should not be assigned to
    REQUIRE(target(x0t-1, y0t) == targetv);
    REQUIRE(target(x0t, y0t-1) == targetv);
    REQUIRE(target(x1t - 1, y0t - 1) == targetv);
    REQUIRE(target(x0t - 1, y0t - 1) == targetv);

    // Test the areas that should be assigned to
    REQUIRE(target(x0t, y0t) != targetv);
    REQUIRE(target(x0t, y0t) == source(x0s, y0s));
    REQUIRE(target(x0t, y1t - 1) == source(x0s, y1s - 1));
    REQUIRE(target(x1t - 1, y0t) == source(x1s - 1, y0s));
    REQUIRE(target(x1t - 1, y1t - 1) == source(x1s - 1, y1s - 1));

    /*****************************************************************/
    // Test negative step
    target = targetv;

    x0s = 0;
    x1s = 3;
    y0s = 3;
    y1s = ny-3;

    x0t = nx - 1;
    x1t = x0t - (x1s - x0s);
    Slice::Int xstep = -1;
    y0t = y1s - 1;
    y1t = y0s - 1;
    Slice::Int ystep = -1;

    ModelArraySlice sourceSlice2(source, {{{x0s, x1s}, {y0s, y1s}}});
    ModelArraySlice targetSlice2(target, {{{x0t, x1t, -1}, {y0t, y1t, -1}}});

    targetSlice2 = sourceSlice2;

    // Test the areas that should not be assigned to
    REQUIRE(target(x1t, y1t+1) == targetv);
    REQUIRE(target(x1t+1, y1t) == targetv);
    REQUIRE(target(x0t, y1t) == targetv);
    REQUIRE(target(x0t, y0t+1) == targetv);

    // Test the areas that should be assigned to
    REQUIRE(target(x0t, y0t) != targetv);
    REQUIRE(target(x0t, y0t) == source(x0s, y0s));
    REQUIRE(target(x0t, y1t - ystep) == source(x0s, y1s - 1));
    REQUIRE(target(x1t - xstep, y0t) == source(x1s - 1, y0s));
    REQUIRE(target(x1t - xstep, y1t - ystep) == source(x1s - 1, y1s - 1));
}

TEST_CASE("ModelArray to Slice")
{
    ModelArray::setDimension(ModelArray::Dimension::X, nx);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny);
    ModelArray::setDimension(ModelArray::Dimension::Z, nz);

    TwoDField source(ModelArray::Type::TWOD);
    source.resize();
    double sourcev = 1.;
    source = sourcev;
    ThreeDField target(ModelArray::Type::THREED);
    target.resize();
    double targetv = -1.;
    target = targetv;

    const auto k = 2;
    // Check the values before the copy
    REQUIRE(target(0, 0, k) == targetv);
    REQUIRE(target(nx-1, 0, k) == targetv);
    REQUIRE(target(0, ny-1, k) == targetv);
    REQUIRE(target(nx-1, ny-1, k) == targetv);
    ModelArraySlice singleLevelSlice(target, {{{{}, {}}, {{}, {}}, {k}}});
    singleLevelSlice = source;

    // Test the areas that should not be assigned to
    REQUIRE(target(0, 0, k-1) == targetv);
    REQUIRE(target(0, 0, k+1) == targetv);
    REQUIRE(target(nx-1, ny-1, k-1) == targetv);
    REQUIRE(target(nx-1, ny-1, k+1) == targetv);

    // Test the areas that should be assigned to
    REQUIRE(target(0, 0, k) == sourcev);
    REQUIRE(target(nx-1, 0, k) == sourcev);
    REQUIRE(target(0, ny-1, k) == sourcev);
    REQUIRE(target(nx-1, ny-1, k) == sourcev);

}

TEST_CASE("Slice to ModelArray")
{
    ModelArray::setDimension(ModelArray::Dimension::X, nx);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny);
    ModelArray::setDimension(ModelArray::Dimension::Z, nz);

    ThreeDField source(ModelArray::Type::THREED);
    source.resize();
    double sourcev = 1.;
    source = sourcev;
    TwoDField target(ModelArray::Type::TWOD);
    target.resize();
    double targetv = -1.;
    target = targetv;

    const auto k = 2;
    // Check the values before the copy
    REQUIRE(target(0, 0) == targetv);
    REQUIRE(target(nx-1, 0) == targetv);
    REQUIRE(target(0, ny-1) == targetv);
    REQUIRE(target(nx-1, ny-1) == targetv);
    ModelArraySlice singleLevelSlice(source, {{{{}, {}}, {{}, {}}, {k}}});
    target = singleLevelSlice;

    // Test the areas that should be assigned to
    REQUIRE(target(0, 0) == sourcev);
    REQUIRE(target(nx-1, 0) == sourcev);
    REQUIRE(target(0, ny-1) == sourcev);
    REQUIRE(target(nx-1, ny-1) == sourcev);

}

TEST_CASE("Index a ModelArray with a Slice")
{
    ModelArray::setDimension(ModelArray::Dimension::X, nx);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny);
    ModelArray::setDimension(ModelArray::Dimension::Z, nz);

    size_t x0 = 4;
    size_t x1 = 9;
    size_t y0 = 6;
    size_t y1 = 14;

    Slice slice {{{x0, x1}, {y0, y1}}};
    TwoDField data(ModelArray::Type::TWOD);
    double v0 = -1.;
    data = v0;
    double v1 = 1.;
    data[slice] = v1;

    // Test the areas that should not be assigned to
    REQUIRE(data(x0 - 1, y0 - 1) == v0);
    REQUIRE(data(x1, y0 - 1) == v0);
    REQUIRE(data(x0 - 1, y1) == v0);
    REQUIRE(data(x1, y1) == v0);
    // Test the areas that should be assigned to
    REQUIRE(data(x0, y0) == v1);
    REQUIRE(data(x0, y1 - 1) == v1);
    REQUIRE(data(x1 - 1, y0) == v1);
    REQUIRE(data(x1 - 1, y1 - 1) == v1);

    // Write the Slice directly into the indexing operator
    x0 = 3;
    x1 = 11;
    y0 = 8;
    y1 = 13;

    v0 = -5.;
    data = v0;
    v1 = 8.;
    data[{{{x0, x1}, {y0, y1}}}] = v1;

    // Test the areas that should not be assigned to
    REQUIRE(data(x0 - 1, y0 - 1) == v0);
    REQUIRE(data(x1, y0 - 1) == v0);
    REQUIRE(data(x0 - 1, y1) == v0);
    REQUIRE(data(x1, y1) == v0);
    // Test the areas that should be assigned to
    REQUIRE(data(x0, y0) == v1);
    REQUIRE(data(x0, y1 - 1) == v1);
    REQUIRE(data(x1 - 1, y0) == v1);
    REQUIRE(data(x1 - 1, y1 - 1) == v1);
}

TEST_CASE("Test buffers")
{
    ModelArray::setDimension(ModelArray::Dimension::X, nx);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny);

    // Test 1 dimensional slices into and out of common buffer types
    // std::vector
    OneDField oned(ModelArray::Type::ONED);
    oned.resize();
    const auto x0 = 3;
    const auto x1 = 11;
    auto oneSlice = oned[{{{x0, x1}}}];

    // Check the functions throw when given too small a buffer
    std::vector<double> shortBuffer;
    shortBuffer.resize(x1-x0-1);
    REQUIRE_THROWS(oneSlice = shortBuffer);
    REQUIRE_THROWS(oneSlice.copyToBuffer(shortBuffer));
    // Fill the array with index values
    for (auto i = 0; i < oned.size(); ++i) {
        oned[i] = i;
    }

    std::vector<double> vectorBuffer;
    vectorBuffer.resize(x1-x0);
    // Fill with index values
    for (auto i = 0; i < vectorBuffer.size(); ++i) {
        vectorBuffer[i] = i;
    }
    // assign from the buffer
    oneSlice = vectorBuffer;
    REQUIRE(oned[x0-1] == x0 - 1);
    REQUIRE(oned[x1] == x1);
    REQUIRE(oned[x0 + 0] == 0);
    REQUIRE(oned[x1 - 1] == x1 - x0 - 1);
    // Refill with index values
    for (auto i = 0; i < oned.size(); ++i) {
        oned[i] = i;
    }
    oneSlice.copyToBuffer(vectorBuffer);
    REQUIRE(vectorBuffer[0] == x0);
    REQUIRE(vectorBuffer[x1 - x0 - 1] == x1 - 1);

    // std::array<x1-x0>
    std::array<double, x1-x0> arrayBuffer;
    for (auto i = 0; i < arrayBuffer.size(); ++i) {
        arrayBuffer[i] = i;
    }
    // assign from the buffer
    oneSlice = arrayBuffer;
    REQUIRE(oned[x0-1] == x0 - 1);
    REQUIRE(oned[x1] == x1);
    REQUIRE(oned[x0 + 0] == 0);
    REQUIRE(oned[x1 - 1] == x1 - x0 - 1);
    // Refill with index values
    for (auto i = 0; i < oned.size(); ++i) {
        oned[i] = i;
    }
    oneSlice.copyToBuffer(arrayBuffer);
    REQUIRE(arrayBuffer[0] == x0);
    REQUIRE(arrayBuffer[x1 - x0 - 1] == x1 - 1);

}

TEST_CASE("General slice creation")
{
    ModelArray::setDimension(ModelArray::Dimension::X, nx);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny);

    size_t x0 = 4;
    size_t x1 = 9;
    size_t y0 = 6;
    size_t y1 = 14;

    TwoDField data(ModelArray::Type::TWOD);
    data.resize();
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            data(i, j) = 100 * j + i;
        }
    }

    // TODO: Make it possible to use a implicit whole-array initalizer here.
    // auto slice = data[{{{},{0}}}];
    auto slice = data[{{{{}}, {0}}}];
    // receiving buffer
    std::vector<double> buff(nx);
    slice.copyToBuffer(buff);
    for (size_t i = 0; i < nx; ++i) {
        REQUIRE(buff[i] == i);
    }
}

TEST_CASE("Iterators")
{
    ModelArray::setDimension(ModelArray::Dimension::X, nx);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny);

    // Equality
    TwoDField a(ModelArray::Type::TWOD);
    ModelArraySlice aSlice(a, {{{3, 7}, {8, 12}}});
    ModelArraySlice::iterator aIter(aSlice);
    ModelArraySlice::iterator aIter2(aSlice);
    REQUIRE(aIter == aIter2);
    ++aIter;
    REQUIRE_FALSE(aIter == aIter2);
    ++aIter2;
    // Begin & end functions
    TwoDField b(ModelArray::Type::TWOD);
    // Fill b
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            b(i, j) = j * 100 + i;
        }
    }

    ModelArraySlice bSlice = b[{{{1, 7}, {9, 11}}}];
    ModelArraySlice::iterator bIter(bSlice);
    REQUIRE(bIter == bSlice.begin());
    REQUIRE(*bIter == 901);
    size_t count = 0;
    while (bIter != bSlice.end()) {
        ++bIter;
        ++count;
    }
    REQUIRE(count == 6 * 2);

    // For loop and non-unit step
    ModelArraySlice bSlice2 = b[{{{2, 9, 2}, {0, 11, 3}}}];
    count = 0;
    for (auto it = bSlice2.begin(); it != bSlice2.end(); ++it) {
        ++count;
    }
    REQUIRE(count == 4 * 4);

    OneDField c(ModelArray::Type::ONED);
    for (size_t i = 0; i < c.size(); ++i) {
        c[i] = i + 100;
    }
    auto cSlice = c[{{{4, 9, 2}}}];
    REQUIRE(*cSlice.begin() == 104);

    for (auto it = cSlice.begin(); it != cSlice.end(); ++it) {
        *it = 0.;
    }
    // Check that earlier values have not been set to 0
    for (size_t i = 0; i < 4; ++i) {
        REQUIRE_FALSE(c[i] == 0);
    }
    // Check the zero values using a for loop over indices
    for (size_t i = 4; i < 9; i += 2) {
        REQUIRE(c[i] == 0.);
    }
    // Check the interleaved values are still non-zero
    for (size_t i = 5; i < 9; i += 2) {
        REQUIRE_FALSE(c[i] == 0);
    }
    // Check the values after the slice are still non-zero
    for (size_t i = 9; i < c.size(); ++i) {
        REQUIRE_FALSE(c[i] == 0);
    }

    count = 0;
    // Check the zero values using the slice iterator
    for (auto& v : cSlice)
    {
        REQUIRE(v == 0.);
        ++count;
    }
    REQUIRE(count == 3);
}

TEST_CASE("Eigen copying")
{
    ModelArray::setDimension(ModelArray::Dimension::X, nx);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny);

    const size_t DG = 6;
    const size_t nWrap = 2;

    Eigen::Matrix<double, Eigen::Dynamic, DG, Eigen::RowMajor> eig;
    eig.resize((nx + 2 * nWrap) * ny, DG);
    SliceIter::MultiDim eigDim = {nx + 2 * nWrap, ny};
    Slice leftColumn {{{0, nWrap}, {}}};
    Slice rightColumn {{{-nWrap, {}}, {}}};
    Slice centreBlock {{{nWrap, -nWrap}, {}}};
    Slice wholeArray{{{}, {}}};
    SliceIter eigLeft(leftColumn, eigDim);
    SliceIter eigRight(rightColumn, eigDim);
    SliceIter eigCentre(centreBlock, eigDim);

    TwoDField source(ModelArray::Type::TWOD);
    source.resize();
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            source(i, j) = i + 100*j;
        }
    }
    TwoDField sink(ModelArray::Type::TWOD);
    sink.resize();
    sink = -1;

    auto eig0 = eig.col(0);
    source[wholeArray].copyToSlicedBuffer(eig0, eigCentre);
    source[leftColumn].copyToSlicedBuffer(eig0, eigRight);
    source[rightColumn].copyToSlicedBuffer(eig0, eigLeft);

    size_t testRow = 3 * ny / 4;
    for (size_t i = 0; i < nx + 2*nWrap; ++i) {
        REQUIRE(eig(Indexer::indexer(eigDim, {i, testRow}), 0) == source ((i + nx - nWrap) % nx, testRow));
    }
    sink[wholeArray].copyFromSlicedBuffer(eig0, eigCentre);
    for (size_t i = 0; i < nx; ++i) {
        REQUIRE(sink(i, testRow) == source(i, testRow));
    }

    /*****************************************************************/
    // Miscellaneous tests
    const size_t oneWrap = 1;
    SliceIter::MultiDim eig1Dim = {nx + 2 * oneWrap, ny + 2 * oneWrap};
    // The source data is the top row of the source array
    Slice topRow1 {{{ }, {ny-1}}};
    Slice topRow2 {{{ }, {-1}}};
    // The data target is the bottom row of the target array, excluding the first and last points in x
    Slice bottomRowBlock1 {{{1, -1}, {0}}};
    SliceIter eig1BottomRowBlock(bottomRowBlock1, eig1Dim);

    // The target array
    Eigen::Matrix<double, Eigen::Dynamic, DG, Eigen::RowMajor> eig1;
    eig1.resize(eig1Dim[0] * eig1Dim[1], DG);
    eig1.setConstant(-1.);

    auto eig1_0 = eig1.col(0);
    source[topRow1].copyToSlicedBuffer(eig1_0, eig1BottomRowBlock);
    // Check the copied values
    size_t iTest = 0;
    REQUIRE(eig1(Indexer::indexer(eig1Dim, {iTest + oneWrap, 0}), 0) == source(iTest, ny-1));

    eig1.setConstant(-1.);

    source[topRow2].copyToSlicedBuffer(eig1_0, eig1BottomRowBlock);
    // Check the copied values
    REQUIRE(eig1(Indexer::indexer(eig1Dim, {iTest + oneWrap, 0}), 0) == source(iTest, ny-1));
}

TEST_CASE("Eigen (ModelArray::DataType) buffers")
{
    ModelArray::setDimension(ModelArray::Dimension::X, nx);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny);

    const size_t DG = 6;
    ModelArray::setDimension(ModelArray::Dimension::COMP, DG);

    size_t sliceNx = 7;
    size_t sliceNy = 11;
    size_t sliceX0 = 3;
    size_t sliceY0 = 5;

    double xMul = 100;
    double yMul = 100;
    ModelArray::DataType eaSrc;
    eaSrc.resize(sliceNx * sliceNy, DG);
    for (size_t j = 0; j < sliceNy; ++j) {
        for (size_t i = 0; i < sliceNx; ++i) {
            size_t ii = Indexer::indexer({sliceNx, sliceNy}, {i, j});
            for (size_t c = 0; c < DG; ++c) {
                eaSrc(ii, c) = c + xMul * (i + yMul * j);
            }
        }
    }
    REQUIRE(eaSrc(Indexer::indexer({sliceNx, sliceNy}, {2, 3}), 4) == 4 + xMul*(2 + yMul * 3));

    ModelArray maInter(ModelArray::Type::TWOCOMP);
    maInter = -1.;
    REQUIRE(maInter.components({sliceX0+2, sliceY0+3})[4] == -1.);
    Slice maSlice {{{sliceX0, sliceX0 + sliceNx}, {sliceY0, sliceY0 + sliceNy}}};
    maInter[maSlice] = eaSrc;
    REQUIRE(maInter.components({sliceX0-1, sliceY0-1})[0] == -1.);
    REQUIRE(maInter.components({sliceX0-1, sliceY0-1})[DG-1] == -1.);
    REQUIRE(maInter.components({sliceX0+0, sliceY0+0})[0] == 0.);
    REQUIRE(maInter.components({sliceX0+0, sliceY0+0})[DG-1] == DG-1);
    REQUIRE(maInter.components({sliceX0+1, sliceY0+0})[0] == xMul);
    REQUIRE(maInter.components({sliceX0+1, sliceY0+0})[DG-1] == xMul + DG-1);
    REQUIRE(maInter.components({sliceX0+0, sliceY0+1})[0] == xMul * yMul);
    REQUIRE(maInter.components({sliceX0+0, sliceY0+1})[DG-1] == xMul  * yMul + DG-1);

    REQUIRE(maInter.components({sliceX0+2, sliceY0+3})[4] == 4 + xMul*(2 + yMul * 3));

    ModelArray::DataType eaSink;
    eaSink.resize(sliceNx * sliceNy, DG);
    eaSink = -1;
    eaSink = maInter[maSlice];
    REQUIRE(eaSink(Indexer::indexer({sliceNx, sliceNy}, {0, 0}), 0) == 0);
    REQUIRE(eaSink(Indexer::indexer({sliceNx, sliceNy}, {0, 0}), 5) == 5);
    REQUIRE(eaSink(Indexer::indexer({sliceNx, sliceNy}, {3, 5}), 2) == 2 + xMul * (3 + yMul * 5));
}

TEST_CASE("Eigen (ModelArray::DataType) subscripted buffers")
{
    ModelArray::setDimension(ModelArray::Dimension::X, nx);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny);

    const size_t DG = 6;
    ModelArray::setDimension(ModelArray::Dimension::COMP, DG);

    double xMul = 100;
    double yMul = 100;

    // Buffer of the array perimeter
    ModelArray maSrc(ModelArray::Type::TWOCOMP);
    maSrc.resize();

    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            double spatial = xMul * (i + yMul * j);
            for (size_t c = 0; c < DG; ++c) {
                maSrc.components({i, j})[c] = spatial + c;
            }
        }
    }
    REQUIRE(maSrc.components({2, 3})[4] == 4 + xMul * (2 + yMul * 3));

    ModelArray::DataType perimeterBuffer;
    perimeterBuffer.resize(2*nx + 2*ny, DG);
    Slice left {{{0}, {}}};
    Slice right {{{nx-1},{}}};
    Slice bottom {{{},{0}}};
    Slice top {{{}, {ny-1}}};
    // bottom 0, nx points
    perimeterBuffer(Eigen::seqN(0, nx), Eigen::all) = static_cast<ModelArray::DataType>(maSrc[bottom]);
    // right nx, ny points
    perimeterBuffer(Eigen::seqN(nx, ny), Eigen::all) = static_cast<ModelArray::DataType>(maSrc[right]);
    // top nx + ny, nx points
    perimeterBuffer(Eigen::seqN(nx+ny, nx), Eigen::all) = static_cast<ModelArray::DataType>(maSrc[top]);
    // left 2nx + ny, ny points
    perimeterBuffer(Eigen::seqN(2*nx+ny, ny), Eigen::all) = static_cast<ModelArray::DataType>(maSrc[left]);

    // bottom check
    for (size_t i = 0; i < nx; ++i) {
        REQUIRE(perimeterBuffer(i, i % DG) == (i % DG) + xMul * i);
    }
    // right check
    for (size_t j = 0; j < ny; ++j) {
        REQUIRE(perimeterBuffer(j + nx, j % DG) == (j % DG) + xMul * (nx-1 + yMul * j));
    }
    // top check
    for (size_t i = 0; i < nx; ++i) {
        REQUIRE(perimeterBuffer(i + nx + ny, i % DG) == (i % DG) + xMul * (i + yMul * (ny -1)));
    }
    // left check
    for (size_t j = 0; j < ny; ++j) {
        REQUIRE(perimeterBuffer(j + 2 * nx + ny, j % DG) == (j % DG) + xMul * (0 + yMul * j));
    }

    ModelArray maSnk(ModelArray::Type::TWOCOMP);
    maSnk.resize();
    maSnk[bottom] = ModelArray::DataType(perimeterBuffer(Eigen::seqN(0, nx), Eigen::all));
    maSnk[right] = ModelArray::DataType(perimeterBuffer(Eigen::seqN(nx, ny), Eigen::all));
    maSnk[top] = ModelArray::DataType(perimeterBuffer(Eigen::seqN(nx + ny, nx), Eigen::all));
    maSnk[left] = ModelArray::DataType(perimeterBuffer(Eigen::seqN(2*nx + ny, ny), Eigen::all));

    // bottom check
    for (size_t i = 0; i < nx; ++i) {
        REQUIRE(maSnk.components({i, 0})[i % DG] == (i % DG) + xMul * i);
    }
    // right check
    for (size_t j = 0; j < ny; ++j) {
        REQUIRE(maSnk.components({nx-1, j})[j % DG] == (j % DG) + xMul * (nx-1 + yMul * j));
    }
    // top check
    for (size_t i = 0; i < nx; ++i) {
        REQUIRE(maSnk.components({i, ny-1})[i % DG] == (i % DG) + xMul * (i + yMul * (ny -1)));
    }
    // left check
    for (size_t j = 0; j < ny; ++j) {
        REQUIRE(maSnk.components({0, j})[j % DG] == (j % DG) + xMul * (0 + yMul * j));
    }

}
TEST_SUITE_END();
} // namespace Nextsim
