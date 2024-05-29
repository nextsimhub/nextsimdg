/*!
 * @file DGModelArray_test.cpp
 *
 * @brief Test that the functions to convert from the dynamics code DGVector
 * to and from ModelArray function correctly.
 *
 * @date Oct 6, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "include/DGModelArray.hpp"

#include "include/ParametricMesh.hpp"

Nextsim::COORDINATES CoordinateSystem = Nextsim::CARTESIAN;

namespace Nextsim {

TEST_SUITE_BEGIN("DGModelArray");
TEST_CASE("DGVector from ModelArray")
{
    static const int DG = 3;
    const size_t nx = 32;
    const size_t ny = 32;
    const size_t mx = 100;
    const size_t my = 100;

    ModelArray::setDimensions(ModelArray::Type::DG, { nx, ny });
    ModelArray::setNComponents(ModelArray::Type::DG, DG);

    ModelArray source(ModelArray::Type::DG);
    source.resize();
    // Fill with data
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; j++) {
            for (size_t c = 0; c < DG; ++c) {
                source.components({i, j})[c] = c + my * (j + mx * (i));
            }
        }
    }

    CHECK(source.components({14, 12})[2] == (2 + my * (12 + mx * (14))));

    // Create the ParametricMesh object
    ParametricMesh smesh(CoordinateSystem);
    smesh.nx = nx;
    smesh.ny = ny;
    smesh.nnodes = nx * ny;
    smesh.nelements = nx * ny;
    smesh.vertices.resize(smesh.nelements, Eigen::NoChange);
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            smesh.vertices(i * ny + j, 0) = i;
            smesh.vertices(i * ny + j, 1) = j;
        }
    }
    DGVector<DG> dest(smesh);
    DGModelArray::ma2dg<DG>(source, dest);

    // Did it work?
    size_t targetPoint = ModelArray::indexFromLocation(ModelArray::Type::DG, {14, 12});
    REQUIRE(dest(targetPoint, 2) == source.components({14, 12})[2]);
}

TEST_CASE("DGVector from ModelArray::Type::H")
{
    static const int DG = 3;
    const size_t nx = 32;
    const size_t ny = 32;
    const size_t mx = 100;
    const size_t my = 100;

    ModelArray::setDimensions(ModelArray::Type::H, { nx, ny });

    ModelArray source(ModelArray::Type::H);
    source.resize();
    // Fill with data
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; j++) {
            source(i, j) = j + mx * (i);
        }
    }

    CHECK(source(14, 12) == 12 + mx * (14));

    // Create the ParametricMesh object
    ParametricMesh smesh(CoordinateSystem);
    smesh.nx = nx;
    smesh.ny = ny;
    smesh.nnodes = nx * ny;
    smesh.nelements = nx * ny;
    smesh.vertices.resize(smesh.nelements, Eigen::NoChange);
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            smesh.vertices(i * ny + j, 0) = i;
            smesh.vertices(i * ny + j, 1) = j;
        }
    }
    DGVector<DG> dest(smesh);
    DGModelArray::ma2dg<DG>(source, dest);

    // Did it work?
    size_t targetPoint = ModelArray::indexFromLocation(ModelArray::Type::H, {14, 12});
    REQUIRE(dest(targetPoint, 0) == source(14, 12));
    REQUIRE(dest(targetPoint, 2) == 0.);

}


TEST_CASE("ModelArray from DGVector")
{
    static const int DG = 6;
    const size_t nx = 31;
    const size_t ny = 33;
    const size_t mx = 100;
    const size_t my = 100;

    // Create the ParametricMesh object
    ParametricMesh smesh(CoordinateSystem);
    smesh.nx = nx;
    smesh.ny = ny;
    smesh.nnodes = nx * ny;
    smesh.nelements = nx * ny;
    smesh.vertices.resize(smesh.nelements, Eigen::NoChange);
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            smesh.vertices(i * ny + j, 0) = i;
            smesh.vertices(i * ny + j, 1) = j;
        }
    }
    DGVector<DG> source(smesh);
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

TEST_CASE("Test with DG = 1") // (It would be a silly case to get wrong!)
{
    static const int DG = 1;
    const size_t nx = 23;
    const size_t ny = 27;
    const size_t mx = 100;

    // Create the ParametricMesh object
    ParametricMesh smesh(CoordinateSystem);
    smesh.nx = nx;
    smesh.ny = ny;
    smesh.nnodes = nx * ny;
    smesh.nelements = nx * ny;
    smesh.vertices.resize(smesh.nelements, Eigen::NoChange);
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            smesh.vertices(i * ny + j, 0) = i;
            smesh.vertices(i * ny + j, 1) = j;
        }
    }
    DGVector<DG> source(smesh);
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

}

// Tests to and from HFields.
TEST_CASE("Test the HField/DG0 transfer") // (It would be a silly case to get wrong!)
{
    static const int DG = 3;
    const size_t nx = 19;
    const size_t ny = 21;
    const size_t mx = 100;
    const size_t my = 100;

    // Create the ParametricMesh object
    ParametricMesh smesh(CoordinateSystem);
    smesh.nx = nx;
    smesh.ny = ny;
    smesh.nnodes = nx * ny;
    smesh.nelements = nx * ny;
    smesh.vertices.resize(smesh.nelements, Eigen::NoChange);
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            smesh.vertices(i * ny + j, 0) = i;
            smesh.vertices(i * ny + j, 1) = j;
        }
    }
    DGVector<DG> source(smesh);
    // Let the H arrays do the index calculation
    ModelArray::setDimensions(ModelArray::Type::H, { nx, ny });
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; j++) {
            size_t index = ModelArray::indexFromLocation(ModelArray::Type::H, {i, j});
            for (size_t c = 0; c < DG; ++c) {
                source(index, c) = j + mx * i + mx * my * c;
            }
        }
    }

    size_t targetPoint = ModelArray::indexFromLocation(ModelArray::Type::H, {14, 12});
    CHECK(source(targetPoint, 0) == 12 + mx * 14);

    ModelArray::setDimensions(ModelArray::Type::DG, { nx, ny });
    ModelArray::setNComponents(ModelArray::Type::DG, DG);

    HField dest(ModelArray::Type::H);
    dest.resize();
    DGModelArray::dg2hField(source, dest);

    // Did it work?
    REQUIRE(dest(14, 12) == source(targetPoint, 0));
    targetPoint = ModelArray::indexFromLocation(ModelArray::Type::H, {12, 14});
    REQUIRE(dest(12, 14) == source(targetPoint, 0));

    // Change the values
    dest *= 2.;
    dest += 1.;

    targetPoint = ModelArray::indexFromLocation(ModelArray::Type::H, {11, 13});
    // Make sure it is different before the copy
    REQUIRE(source(targetPoint, 0) != dest(11, 13));
    double retained = source(targetPoint, 1);
    REQUIRE(source(targetPoint, 0) != retained);

    DGModelArray::hField2dg<DG>(dest, source);
    // And identical again after
    REQUIRE(source(targetPoint, 0) == dest(11, 13));
    REQUIRE(source(targetPoint, 1) == retained);
    REQUIRE(source(targetPoint, 0) != retained);

}
TEST_SUITE_END();

}
