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

TEST_CASE("CGVector from ModelArray", "[DGModelArray]")
{
    static const int CG = 2;
    // Base grid size
    const size_t nx = 37;
    const size_t ny = 31;
    const size_t mx = 100;
    const size_t my = 100;

    ModelArray::setDimensions(ModelArray::Type::CG, { CG * nx + 1, CG * ny + 1 });

    ModelArray source(ModelArray::Type::CG);
    source.resize();
    // Fill with data
    for (size_t i = 0; i < CG * nx + 1; ++i) {
        for (size_t j = 0; j < CG * ny + 1; j++) {
                source(i, j) = j + mx * i;
            }
        }

    CHECK(source(54, 42) == (54 * mx + 42));

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
    CGVector<CG> dest(smash);
    CGModelArray::ma2cg<CG>(source, dest);

    // Did it work?
    size_t targetPoint = ModelArray::indexFromLocation(ModelArray::Type::CG, {52, 44});
    REQUIRE(dest(targetPoint) == source(52, 44));
}
/*
TEST_CASE("ModelArray from DGVector", "[DGModelArray]")
{
    static const int DG = 6;
    const size_t nx = 31;
    const size_t ny = 33;
    const size_t mx = 100;
    const size_t my = 100;

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
    DGVector<DG> source(smash);
    // Let the H arrays do the index calculation
    ModelArray::setDimensions(ModelArray::Type::H, { nx, ny });
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; j++) {
            size_t index = ModelArray::indexFromLocation(ModelArray::Type::H, {i, j});
            for (size_t c = 0; c < DG; ++c) {
                source(index, c) = c + my * (j + mx * (i));
            }
        }
    }
    size_t targetPoint = ModelArray::indexFromLocation(ModelArray::Type::H, {14, 12});
    CHECK(source(targetPoint, 4) == 4 + my * (12 + mx *(14)));

    ModelArray::setDimensions(ModelArray::Type::DG, { nx, ny });
    ModelArray::setNComponents(ModelArray::Type::DG, DG);

    ModelArray dest(ModelArray::Type::DG);
    dest.resize();
    DGModelArray::dg2ma(source, dest);

    // Did it work?
    REQUIRE(dest.components({14, 12})[4] == source(targetPoint, 4));
}

TEST_CASE("Test with DG = 1", "[DGModelArray]") // (It would be a silly case to get wrong!)
{
    static const int DG = 1;
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
    DGVector<DG> source(smash);
    // Let the H arrays do the index calculation
    ModelArray::setDimensions(ModelArray::Type::H, { nx, ny });
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; j++) {
            size_t index = ModelArray::indexFromLocation(ModelArray::Type::H, {i, j});
            for (size_t c = 0; c < DG; ++c) {
                source(index, c) = j + mx * i;
            }
        }
    }

    size_t targetPoint = ModelArray::indexFromLocation(ModelArray::Type::H, {14, 12});
    CHECK(source(targetPoint) == 12 + mx * 14);

    ModelArray::setDimensions(ModelArray::Type::DG, { nx, ny });
    ModelArray::setNComponents(ModelArray::Type::DG, DG);

    ModelArray dest(ModelArray::Type::DG);
    dest.resize();
    DGModelArray::dg2ma(source, dest);

    // Did it work?
    REQUIRE(dest(14, 12) == source(targetPoint));
    targetPoint = ModelArray::indexFromLocation(ModelArray::Type::H, {12, 14});
    REQUIRE(dest(12, 14) == source(targetPoint));

    // Change the values
    dest *= 2.;
    dest += 1.;

    targetPoint = ModelArray::indexFromLocation(ModelArray::Type::H, {11, 13});
    // Make sure it is different before the copy
    REQUIRE(source(targetPoint) != dest(11, 13));
    DGModelArray::ma2dg<DG>(dest, source);
    // And identical again after
    REQUIRE(source(targetPoint) == dest(11, 13));

}*/
}
