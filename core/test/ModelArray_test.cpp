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
    // Check neighbouring y indices differ in value by 1
    REQUIRE(check1d(x, y) - check1d(x, y-1) == 1);
    // Check neighbouring x values differ in value by ny
    REQUIRE(check1d(x, y) - check1d(x-1, y) == dims2[1]);

    REQUIRE(check1d(dims2[0]-1, dims2[1]-1) == dims2[0] * dims2[1] - 1);
}

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
    REQUIRE(check4d(4, 7, 2, 5) == 4725);

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

    REQUIRE(check4d(4, 7, 2, 5) == 4725);

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

    REQUIRE(check4d(4, 7, 2, 5) == 4725);


    REQUIRE(check4d[{4, 7, 2, 6}] == 4726);
}

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

    size_t target = (((i) * dims4[1] + j) * dims4[2] + k) * dims4[3] + l;
    REQUIRE(primorial[target] == target);

    REQUIRE(primorial(i, j, k, l) == target);

}

TEST_CASE("Moving data")
{
    size_t n = 10;
    ModelArray::setDimensions(ModelArray::Type::TWOD, {n, n});

    ModelArray src = ModelArray::TwoDField();
    for (int i = 0; i < n * n; ++i) {
        src[i] = i;
    }

    ModelArray cpyCtor(src);
    REQUIRE(cpyCtor(2, 3) == src(2, 3));

    ModelArray cpyAss = ModelArray::TwoDField();
    cpyAss = src;
    REQUIRE(cpyAss(2, 3) == 23);
}

TEST_CASE("Instance setDimensions sets instance dimensions")
{
    DosDField uu = ModelArray::DosDField();
    ModelArray::MultiDim udim = {5, 5};
    uu.setDimensions(udim);
    REQUIRE(uu.size() == udim[0] * udim[1]);
    REQUIRE(uu.nDimensions() == 2);
    REQUIRE(uu.dimensions() == udim);
}

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
TEST_SUITE_END();

} /* namespace Nextsim */
