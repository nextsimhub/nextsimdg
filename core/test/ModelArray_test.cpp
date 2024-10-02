/*!
 * @file ModelData_test.cpp
 *
 * @date Feb 24, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "include/ModelArray.hpp"

namespace Nextsim {

TEST_SUITE_BEGIN("ModelArray");
// Test that the (special case) two dimensional index functions correctly
TEST_CASE("Two dimensional data access test")
{
    ModelArray::MultiDim dims2 = {15, 25};

    ModelArray::setDimensions(ModelArray::Type::TWOD, dims2);

    ModelArray check1d = ModelArray::TwoDField();

    REQUIRE(check1d.nDimensions() == 2);

    // Check the ordering: later indexes should vary faster
    // Fill data as 1d
    for (size_t i = 0; i < dims2[0] * dims2[1]; ++i) {
        check1d[i] = i;
    }

    size_t x = 7;
    size_t y = 13;
    // Check neighbouring x indices differ in value by 1
    REQUIRE(check1d(x+1, y) - check1d(x, y) == 1);
    // Check neighbouring y indices differ in value by nx
    REQUIRE(check1d(x, y+1) - check1d(x, y) == dims2[0]);

    REQUIRE(check1d(dims2[0]-1, dims2[1]-1) == dims2[0] * dims2[1] - 1);
}

// Test that higher dimensional indexing functions correctly
TEST_CASE("Higher dimensional indexing")
{
    size_t dimLen = 10;
    size_t arrayLen = dimLen * dimLen * dimLen * dimLen;
    ModelArray::MultiDim dims4 = {dimLen, dimLen, dimLen, dimLen};
    ModelArray::setDimensions(ModelArray::Type::FOURD, dims4);

    ModelArray check4d = ModelArray::FourDField();

    REQUIRE(check4d.nDimensions() == 4);
    REQUIRE(check4d.size() == dimLen * dimLen * dimLen * dimLen);

    for (size_t i = 0; i < check4d.size(); ++i) {
        check4d[i] = i;
    }

    // Check indexing using the fact that the dimensions are the same as our counting base
    REQUIRE(check4d(4, 7, 2, 5) == 5274);

    // Reset the data to zero
    for (size_t i = 0; i < check4d.size(); ++i) {
        check4d[i] = 0.;
    }
    REQUIRE(check4d(4, 7, 2, 5) == 0);

    double* data = new double[arrayLen];
    for (size_t i = 0; i < arrayLen; ++i) {
        data[i] = i;
    }

    check4d.setData(data);

    REQUIRE(check4d(4, 7, 2, 5) == 5274);

    // Reset the data to zero
    for (size_t i = 0; i < check4d.size(); ++i) {
        check4d[i] = 0.;
    }
    REQUIRE(check4d(4, 7, 2, 5) == 0);

    std::vector<double> vData = std::vector<double>(arrayLen);
    for (size_t i = 0; i < arrayLen; ++i) {
        vData[i] = i;
    }

    check4d.setData(vData.data());

    REQUIRE(check4d(4, 7, 2, 5) == 5274);


    REQUIRE(check4d[{5, 7, 2, 5}] == 5275);
}

// Test that higher dimensional indexing functions correctly
TEST_CASE("Higher dimensional indexing 2")
{
    ModelArray::MultiDim dims4 = {3, 5, 7, 11};
    size_t totalSize = dims4[0] * dims4[1] * dims4[2] * dims4[3];
    ModelArray::setDimensions(ModelArray::Type::FOURD, dims4);

    FourDField primorial = ModelArray::FourDField();

    REQUIRE(primorial.nDimensions() == 4);
    REQUIRE(primorial.size() == totalSize);

    for (size_t i = 0; i < primorial.size(); ++i) {
        primorial[i] = i;
    }

    size_t i = 2;
    size_t j = 4;
    size_t k = 5;
    size_t l = 7;

    size_t target = i + dims4[0] * (j + dims4[1] * (k + dims4[2] * (l)));
            //(((i) * dims4[1] + j) * dims4[2] + k) * dims4[3] + l;
    REQUIRE(primorial[target] == target);

    REQUIRE(primorial(i, j, k, l) == target);

}

// Test that the copy constructor and copy assignment operator initialize that
// data correctly.
TEST_CASE("Copy constructor and copy assignment operator")
{
    size_t n = 10;
    ModelArray::setDimensions(ModelArray::Type::TWOD, {n, n});

    ModelArray src = ModelArray::TwoDField();
    for (int i = 0; i < n * n; ++i) {
        src[i] = i;
    }

    // Test the copy constructor
    ModelArray copyConstructor(src);
    REQUIRE(copyConstructor(2, 3) == src(2, 3));

    // Test copy assignment
    ModelArray copyAssignment = ModelArray::TwoDField();
    copyAssignment = src;
    REQUIRE(copyAssignment(2, 3) == src(2, 3));
}

// Test that setting the dimension via the function applied to an instance
// correctly propagates to the dimensions of the type.
TEST_CASE("Instance setDimensions sets instance dimensions")
{
    ZUField uu = ModelArray::ZUField();
    ModelArray::MultiDim udim = {5, 5};
    uu.setDimensions(udim);
    REQUIRE(uu.size() == udim[0] * udim[1]);
    REQUIRE(uu.nDimensions() == 2);
    REQUIRE(uu.dimensions() == udim);
}

// Test the arithmetic operators of the class.
TEST_CASE("Arithmetic tests")
{
    // Only test HField for now
    ModelArray::setDimensions(ModelArray::Type::ONED, {2});
    OneDField lhs;
    OneDField rhs;
    lhs[0] = 9.;
    lhs[1] = 10.;
    rhs[0] = 3.;
    rhs[1] = -5.;

    OneDField sum = lhs + rhs;
    REQUIRE(sum[0] == 12.);
    REQUIRE(sum[1] == 5.);
    OneDField difference = lhs - rhs;
    REQUIRE(difference[0] == 6.);
    REQUIRE(difference[1] == 15.);
    OneDField product = lhs * rhs;
    REQUIRE(product[0] == 27.);
    REQUIRE(product[1] == -50.);
    OneDField quotient = lhs / rhs;
    REQUIRE(quotient[0] == 3.);
    REQUIRE(quotient[1] == -2.);
    OneDField negative = -rhs;
    REQUIRE(negative[0] == -3);
    REQUIRE(negative[1] == 5);

    double three = 3;
    double four = 4;
    sum = lhs + three;
    REQUIRE(sum[0] == 12);
    REQUIRE(sum[1] == 13);
    sum = four + rhs;
    REQUIRE(sum[0] == 7);
    REQUIRE(sum[1] == -1);
    difference = lhs - three;
    REQUIRE(difference[0] == 6);
    REQUIRE(difference[1] == 7);
    difference = four - rhs;
    REQUIRE(difference[0] == 1);
    REQUIRE(difference[1] == 9);
    product = lhs * three;
    REQUIRE(product[0] == 27);
    REQUIRE(product[1] == 30);
    product = four * rhs;
    REQUIRE(product[0] == 12);
    REQUIRE(product[1] == -20);
    quotient = lhs / three;
    REQUIRE(quotient[0] == 3);
    REQUIRE(quotient[1] == (10. / 3.));
    quotient = four / rhs;
    REQUIRE(quotient[0] == (4. / 3.));
    REQUIRE(quotient[1] == (4. / -5.));

    // Finally, in-place arithmetic
    lhs += rhs;
    REQUIRE(lhs[0] == 12);
    REQUIRE(lhs[1] == 5);
    lhs -= rhs;
    REQUIRE(lhs[0] == 9);
    REQUIRE(lhs[1] == 10);
    lhs *= rhs;
    REQUIRE(lhs[0] == 27);
    REQUIRE(lhs[1] == -50);
    lhs /= rhs;
    REQUIRE(lhs[0] == 9);
    REQUIRE(lhs[1] == 10);

    lhs += three;
    REQUIRE(lhs[0] == 12);
    REQUIRE(lhs[1] == 13);
    lhs -= four;
    REQUIRE(lhs[0] == 8);
    REQUIRE(lhs[1] == 9);
    lhs *= three;
    REQUIRE(lhs[0] == 24);
    REQUIRE(lhs[1] == 27);
    lhs /= four;
    REQUIRE(lhs[0] == 6);
    REQUIRE(lhs[1] == 6.75);

    OneDField fill = ModelArray::OneDField();
    double filldub = 5.2354;
    fill = filldub;

    REQUIRE(fill[0] == filldub);
    REQUIRE(fill[1] == filldub);
}

// Location from index. Index from location is assumed to work as it is a
// wrapper around indexr()
TEST_CASE("Location from index")
{
    const size_t nx = 31;
    const size_t ny = 37;
    const size_t nz = 41;

    ModelArray::setDimensions(ModelArray::Type::FOURD, {nx, ny, nz, 1});
    size_t x = 13;
    size_t y = 17;
    size_t z = 19;

    size_t index = ModelArray::indexFromLocation(ModelArray::Type::FOURD, {x, y, z, 0});
    ModelArray::MultiDim loc = ModelArray::locationFromIndex(ModelArray::Type::FOURD, index);
    REQUIRE(loc[0] == x);
    REQUIRE(loc[1] == y);
    REQUIRE(loc[2] == z);
}

// Test the zIndexAndLayer function to ensure that it accesses the correct
// point in a three-dimensional ModelArray.
TEST_CASE("zIndexAndLayer")
{
    const size_t nx = 29;
    const size_t ny = 23;
    const size_t nz = 11;

    ModelArray::setDimensions(ModelArray::Type::THREED, {nx, ny, nz});

    ThreeDField threeD(ModelArray::Type::THREED);
    threeD.resize();

    size_t mul = 100;
    // Fill the array fastest last
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t k = 0; k < nz; ++k) {
                threeD(i, j, k) = k + mul * (j + mul * (i));
            }
        }
    }

    size_t x = 19;
    size_t y = 17;
    size_t z = 7;
    size_t ind = ModelArray::indexFromLocation(ModelArray::Type::TWOD, {x, y});
    REQUIRE(threeD.zIndexAndLayer(ind, z) == threeD(x, y, z));

}

// Does zeroing a 3D array actually set all elements to zero?
TEST_CASE("Zero a ThreeDField")
{
    const size_t nx = 23;
    const size_t ny = 19;
    const size_t nz = 5;

    ModelArray::setDimensions(ModelArray::Type::THREED, {nx, ny, nz});

    REQUIRE(ModelArray::size(ModelArray::Dimension::X) == nx);
    REQUIRE(ModelArray::size(ModelArray::Dimension::Y) == ny);
    REQUIRE(ModelArray::size(ModelArray::Dimension::Z) == nz);

    ThreeDField treed(ModelArray::Type::THREED);
    treed.resize();

    auto dims = treed.dimensions();
    REQUIRE(dims[0] == nx);
    REQUIRE(dims[1] == ny);
    REQUIRE(dims[2] == nz);

    REQUIRE(treed.trueSize() == nx * ny * nz);

    treed = 0.;

    double absSum = 0.;
    // Sum the absolute value of each element. This is only equal to zero if every element is equal to zero.
    for (size_t i = 0; i < treed.trueSize(); ++i) {
        absSum += std::fabs(treed[i]);
    }
    REQUIRE(absSum == 0.);
}

// Copy a block of contiguous data, specifically a 2d slice of a 3d array
TEST_CASE("Copying contiguous data")
{
    const size_t nx = 17;
    const size_t ny = 23;
    const size_t nz = 3;

    ModelArray::setDimensions(ModelArray::Type::THREED, {nx, ny, nz});

    // A 3D array and a 2D array
    TwoDField twod(ModelArray::Type::TWOD);
    ThreeDField threed(ModelArray::Type::THREED);

    // Thoroughly check the dimensions
    REQUIRE(ModelArray::size(ModelArray::Dimension::X) == nx);
    REQUIRE(ModelArray::size(ModelArray::Dimension::Y) == ny);
    REQUIRE(ModelArray::size(ModelArray::Dimension::Z) == nz);

    REQUIRE(twod.dimensions()[0] == nx);
    REQUIRE(twod.dimensions()[1] == ny);

    REQUIRE(threed.dimensions()[0] == nx);
    REQUIRE(threed.dimensions()[1] == ny);
    REQUIRE(threed.dimensions()[2] == nz);

    twod.resize();
    threed.resize();
    threed = 0.;

    // Check the value of every element is zero
    double absSum;
    absSum = 0;
    for (size_t idx = 0; idx < threed.size(); ++idx) {
        absSum += std::fabs(threed[idx]);
    }
    REQUIRE(absSum == 0.);

    double dx = 100;

    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            twod(i, j) = dx * i + j;
        }
    }

    // Assign data to the middle z level
    ModelArray::MultiDim dims = threed.dimensions();
    size_t size2D = dims[0] * dims[1];
    threed.setData(twod.data(), 1 * size2D, size2D);

    // Check the first level is all zeros
    double lvl0Sum = 0;
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            lvl0Sum += std::fabs(threed(i, j, 0UL));
        }
    }
    REQUIRE(lvl0Sum == 0.);

    // Check the first level matches the 2D array
    double lvl1Sum = 0;
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            lvl1Sum += std::fabs(threed(i, j, 1UL) - twod(i, j));
        }
    }
    REQUIRE(lvl1Sum == 0.);

    // Check the third level is all zeros
    double lvl2Sum = 0;
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            lvl2Sum += std::fabs(threed(i, j, 2UL));
        }
    }
    REQUIRE(lvl2Sum == 0.);


}
// Test the setData function with source and target start indices
TEST_CASE("Copying contiguous data II")
{
    const size_t nx = 13;
    const size_t ny = 19;
    const size_t nz = 7;

    ModelArray::setDimensions(ModelArray::Type::THREED, {nx, ny, nz});

    // A 3D array and a 2D array
    TwoDField twod(ModelArray::Type::TWOD);
    ThreeDField threed(ModelArray::Type::THREED);

    // Thoroughly check the dimensions
    REQUIRE(ModelArray::size(ModelArray::Dimension::X) == nx);
    REQUIRE(ModelArray::size(ModelArray::Dimension::Y) == ny);
    REQUIRE(ModelArray::size(ModelArray::Dimension::Z) == nz);

    REQUIRE(twod.dimensions()[0] == nx);
    REQUIRE(twod.dimensions()[1] == ny);

    REQUIRE(threed.dimensions()[0] == nx);
    REQUIRE(threed.dimensions()[1] == ny);
    REQUIRE(threed.dimensions()[2] == nz);

    twod.resize();
    threed.resize();

    double twod0 = 2.;
    double threed0 = 3.;

    twod = twod0;
    threed = threed0;

    // Shorter source array, filling z-level 3
    size_t size2D = twod.size();
    size_t k = 3;
    threed.setData(twod, 0, 3 * size2D);

    REQUIRE(threed(0UL, 0UL, 0UL) == threed0);
    REQUIRE(threed(nx - 1, ny - 1, k - 1) == threed0);
    REQUIRE(threed(0UL, 0UL, k) == twod0);
    REQUIRE(threed(nx / 2, ny / 2, k) == twod0);
    REQUIRE(threed(nx - 1, ny - 1, k) == twod0);
    REQUIRE(threed(0UL, 0UL, k + 1) == threed0);
    REQUIRE(threed(nx - 1, ny - 1, nz - 1) == threed0);

    // Shorter target array, filling from z-level 3
    twod0 = 4.;
    threed0 = 5.;
    twod = twod0;

    for (size_t kk = 0; kk < nz; ++kk) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {
                threed(i, j, kk) = threed0 + kk;
            }
        }
    }
    REQUIRE(twod(0UL, 0UL) == twod0);
    REQUIRE(threed(0UL, 0UL, 4UL) == threed0 + 4);
    twod.setData(threed, 3 * size2D, 0);
    REQUIRE(twod(0UL, 0UL) == threed0 + k);
    REQUIRE(twod(0UL, 0UL) == threed0 + k);

    // Fill something that's not a z-level
    twod0 = 6.;
    threed0 = 7.;
    twod = twod0;
    threed = threed0;

    size_t sourceStart = size2D / 2;
    size_t targetStart = threed.indexFromLocation({11, 13, 5});
    size_t targetEnd = targetStart + size2D - sourceStart;
    threed.setData(twod, sourceStart, targetStart);
    REQUIRE(threed[targetStart - 1] == threed0);
    REQUIRE(threed[targetStart] == twod0);
    REQUIRE(threed[targetEnd - 1] == twod0);
    REQUIRE(threed[targetEnd] == threed0);

}
TEST_SUITE_END();

} /* namespace Nextsim */
