/*!
 * @file ParametricMesh_test.cpp
 *
 * @brief Test the ParametricMesh class, especially processing from ModelArray files.
 *
 * @date Dec 15, 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "include/ParametricMesh.hpp"

#include "../test/FakeSmeshData.hpp"
#include "include/gridNames.hpp"

#include <filesystem>

namespace Nextsim {

COORDINATES CoordinateSystem = CARTESIAN;

#define TO_STR(s) TO_STRI(s)
#define TO_STRI(s) #s
#ifndef TEST_FILE_SOURCE
#define TEST_FILE_SOURCE .
#endif

static const std::string smeshFile = TO_STR(TEST_FILE_SOURCE) + std::string("/") + "25km_NH.smesh";
static const size_t nx = 154;
static const size_t ny = 121;
static const double dx = 25000;
static const double dy = 25000;

TEST_SUITE_BEGIN("ParametricMesh");
// Read data from a .smesh file. This is the old way, but also provides a
// reference for testing the newer ways.
TEST_CASE("Test readmesh")
{
    ParametricMesh smesh(CoordinateSystem);

    smesh.readmesh(smeshFile);

    // Sizes of things
    REQUIRE(smesh.nx == nx);
    REQUIRE(smesh.ny == ny);
    REQUIRE(smesh.nelements == nx * ny);
    REQUIRE(smesh.nnodes == (nx + 1) * (ny + 1));
    REQUIRE(smesh.vertices.col(0).size() == (nx + 1) * (ny + 1));
    // Coordinate values
    REQUIRE(smesh.vertices(0, 0) == 0.);
    REQUIRE(smesh.vertices(0, 1) == 0.);
    REQUIRE(smesh.vertices(1, 0) == dx);
    REQUIRE(smesh.vertices(1, 1) == 0.);
    REQUIRE(smesh.vertices(nx + 1, 0) == 0.);
    REQUIRE(smesh.vertices(nx + 1, 1) == dy);
    REQUIRE(smesh.vertices((nx + 1) * (ny + 1) - 1, 0) == nx * dx);
    REQUIRE(smesh.vertices((nx + 1) * (ny + 1) - 1, 1) == ny * dy);

    // Landmask values
    REQUIRE(!smesh.landmask[0]);
    REQUIRE(!smesh.landmask[4]);
    REQUIRE(smesh.landmask[5]);
    REQUIRE(smesh.landmask[nx - 1]);
    REQUIRE(!smesh.landmask[(ny - 1) * nx]);

    // No periodic boundary conditions
    REQUIRE(smesh.periodic.size() == 0);
}

// Read data from ModelArrays and compare to that obtained from a .smesh file.
// The class FakeSmeshData provide ModelArrays containing data equivalent to
// that obtained from read the reference smesh file.
TEST_CASE("Compare readmesh and landmask reading")
{

    ParametricMesh fromFile(CoordinateSystem);
    ParametricMesh fromArrays(CoordinateSystem);

    fromFile.readmesh(smeshFile);

    ModelState fakeSmeshData = FakeSmeshData::getData();
    REQUIRE(fakeSmeshData.data.at(xName).trueSize() == nx * ny);
    REQUIRE(fakeSmeshData.data.at(coordsName).trueSize() == (nx + 1) * (ny + 1));
    REQUIRE(fakeSmeshData.data.at(maskName).trueSize() == nx * ny);

    fromArrays.coordinatesFromModelArray(fakeSmeshData.data.at(coordsName));
    // Sizes of things
    REQUIRE(fromArrays.nx == fromFile.nx);
    REQUIRE(fromArrays.ny == fromFile.ny);
    REQUIRE(fromArrays.nelements == fromFile.nelements);
    REQUIRE(fromArrays.nnodes == fromFile.nnodes);
    REQUIRE(fromArrays.vertices.col(0).size() == fromFile.vertices.col(0).size());
    // Coordinate values
    std::vector<size_t> checkIndices = { 0, 1, nx + 1, (nx + 1) * (ny + 1) - 1 };
    for (auto index : checkIndices) {
        REQUIRE(fromArrays.vertices(index, 0) == fromFile.vertices(index, 0));
        REQUIRE(fromArrays.vertices(index, 1) == fromFile.vertices(index, 1));
    }

    // Landmask values
    fromArrays.landmaskFromModelArray(fakeSmeshData.data.at(maskName));
    for (size_t idx = 0; idx < fromArrays.nelements; ++idx) {
        REQUIRE(fromArrays.landmask[idx] == fromFile.landmask[idx]);
    }

    // No periodic boundary conditions
    REQUIRE(fromArrays.periodic.size() == 0);
}
}
