/*!
 * @file DGModelArray_test.cpp
 *
 * @date Oct 6, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "include/CGModelArray.hpp"

#include "include/ParametricMesh.hpp"

namespace Nextsim {

TEST_CASE("cgDims test", "[CGModelArray]")
{
    static const int CG = 2;
    ModelArray::Dimensions hDims = { 23, 29 };
    ModelArray::setDimensions(ModelArray::Type::H, hDims);
    ModelArray::Dimensions cgDims = CGModelArray::cgDimensions<CG>(hDims);
    REQUIRE(cgDims[0] == CG * hDims[0] + 1);
    REQUIRE(cgDims[1] == CG * hDims[1] + 1);
}

TEST_CASE("CGVector from ModelArray", "[DGModelArray]")
{
    static const int CG = 2;
    // Base grid size
    const size_t nx = 37;
    const size_t ny = 31;
    const size_t mx = 100;

    ModelArray::setDimensions(ModelArray::Type::CG, CGModelArray::cgDimensions<CG>({ nx, ny }));

    ModelArray maSource(ModelArray::Type::CG);
    maSource.resize();
    // Fill with data
    for (size_t i = 0; i < CG * nx + 1; ++i) {
        for (size_t j = 0; j < CG * ny + 1; j++) {
            maSource(i, j) = j + mx * i;
        }
    }

    CHECK(maSource(54, 42) == (54 * mx + 42));

    // Create the ParametricMesh object
    ParametricMesh smash;
    smash.nx = nx;
    smash.ny = ny;
    smash.nnodes = nx * ny;
    smash.nelements = nx * ny;
    smash.vertices.resize(smash.nelements, Eigen::NoChange);
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            smash.vertices(i * ny + j, 0) = i;
            smash.vertices(i * ny + j, 1) = j;
        }
    }
    CGVector<CG> cgDest(smash);
    CGModelArray::ma2cg<CG>(maSource, cgDest);

    // Did it work?
    size_t targetPoint = ModelArray::indexFromLocation(ModelArray::Type::CG, { 52, 44 });
    REQUIRE(cgDest(targetPoint) == maSource(52, 44));
}

TEST_CASE("ModelArray from CGVector", "[CGModelArray]")
{
    static const int CG = 2;
    const size_t nx = 31;
    const size_t ny = 33;
    const size_t mx = 100;

    // Create the ParametricMesh object
    ParametricMesh smash;
    smash.nx = nx;
    smash.ny = ny;
    smash.nnodes = nx * ny;
    smash.nelements = nx * ny;
    smash.vertices.resize(smash.nelements, Eigen::NoChange);
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            smash.vertices(i * ny + j, 0) = i;
            smash.vertices(i * ny + j, 1) = j;
        }
    }
    CGVector<CG> cgSource(smash);
    ModelArray::setDimensions(ModelArray::Type::CG, CGModelArray::cgDimensions<CG>({ nx, ny }));
    for (size_t i = 0; i < CG * nx + 1; ++i) {
        for (size_t j = 0; j < CG * ny + 1; j++) {
            size_t index = ModelArray::indexFromLocation(ModelArray::Type::CG, { i, j });
            cgSource(index) = mx * i + j;
        }
    }
    size_t targetPoint = ModelArray::indexFromLocation(ModelArray::Type::CG, { 14, 12 });
    CHECK(cgSource(targetPoint) == 14 * mx + 12);

    ModelArray maDest(ModelArray::Type::CG);
    maDest.resize();
    CGModelArray::cg2ma(cgSource, maDest);

    // Did it work?
    REQUIRE(maDest(14, 12) == cgSource(targetPoint));
}

TEST_CASE("Test with CG = 1", "[DGModelArray]") // (It would be a silly case to get wrong!)
{
    static const int CG = 1;
    const size_t nx = 23;
    const size_t ny = 27;
    const size_t mx = 100;

    // Create the ParametricMesh object
    ParametricMesh smash;
    smash.nx = nx;
    smash.ny = ny;
    smash.nnodes = nx * ny;
    smash.nelements = nx * ny;
    smash.vertices.resize(smash.nelements, Eigen::NoChange);
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            smash.vertices(i * ny + j, 0) = i;
            smash.vertices(i * ny + j, 1) = j;
        }
    }
    CGVector<CG> cgSource(smash);
    // Let the H arrays do the index calculation
    ModelArray::setDimensions(ModelArray::Type::H, { nx, ny });
    ModelArray::setDimensions(ModelArray::Type::CG,
        CGModelArray::cgDimensions<CG>(ModelArray::dimensions(ModelArray::Type::H)));
    for (size_t i = 0; i < CG * nx + 1; ++i) {
        for (size_t j = 0; j < CG * ny + 1; j++) {
            size_t index = ModelArray::indexFromLocation(ModelArray::Type::CG, { i, j });
            cgSource(index) = mx * i + j;
        }
    }

    size_t targetPoint = ModelArray::indexFromLocation(ModelArray::Type::CG, { 14, 22 });
    CHECK(cgSource(targetPoint) == mx * 14 + 22);

    ModelArray maDest(ModelArray::Type::CG);
    maDest.resize();
    CGModelArray::cg2ma(cgSource, maDest);

    // Did it work?
    REQUIRE(maDest(14, 22) == cgSource(targetPoint));
    targetPoint = ModelArray::indexFromLocation(ModelArray::Type::CG, { 22, 14 });
    REQUIRE(maDest(22, 14) == cgSource(targetPoint));

    // Change the values
    maDest *= 2.;
    maDest += 1.;

    targetPoint = ModelArray::indexFromLocation(ModelArray::Type::CG, { 21, 13 });
    // Make sure it is different before the copy
    REQUIRE(cgSource(targetPoint) != maDest(21, 13));
    CGModelArray::ma2cg<CG>(maDest, cgSource);
    // And identical again after
    REQUIRE(cgSource(targetPoint) == maDest(21, 13));
}
}
