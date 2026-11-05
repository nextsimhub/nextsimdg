/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @brief Test the ParametricMesh class, especially processing from ModelArray files.
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
static const FloatType dx = 25000;
static const FloatType dy = 25000;

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

    // Cached element areas must match a direct recomputation from the vertices
    REQUIRE(smesh.elementArea.size() == smesh.nelements);
    for (size_t i = 0; i < smesh.nelements; ++i) {
        REQUIRE(smesh.elementArea[i] == smesh.compute_area(i));
        REQUIRE(smesh.area(i) == smesh.elementArea[i]);
    }
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

    // Cached element areas must be populated and identical to those from readmesh
    REQUIRE(fromArrays.elementArea.size() == fromFile.nelements);
    for (size_t idx = 0; idx < fromArrays.nelements; ++idx) {
        REQUIRE(fromArrays.area(idx) == fromFile.area(idx));
        REQUIRE(fromArrays.elementArea[idx] == fromArrays.compute_area(idx));
    }
}

// Check that the area cache is correctly recomputed after the geometry
// transformations used in the production initialisation path
// (DynamicsKernel::initialise: coordinatesFromModelArray + RotatePoleToGreenland).
TEST_CASE("Area cache updated by mesh transformations")
{
    ParametricMesh smesh(SPHERICAL);

    // Spherical grid in radians covering the whole globe
    smesh.nx = nx;
    smesh.ny = ny;
    smesh.nelements = nx * ny;
    smesh.nnodes = (nx + 1) * (ny + 1);
    smesh.vertices.resize(smesh.nnodes, 2);
    for (size_t iy = 0; iy <= ny; ++iy)
        for (size_t ix = 0; ix <= nx; ++ix) {
            smesh.vertices((nx + 1) * iy + ix, 0) = -1.0_ft * pi + 2.0_ft * pi * ix / nx;
            smesh.vertices((nx + 1) * iy + ix, 1) = -0.5_ft * pi + 1.0_ft * pi * iy / ny;
        }
    smesh.landmask.resize(nx * ny, true);
    smesh.precalculate_area();

    // Cache must match a direct recomputation, and total area must match the
    // surface of the earth up to machine accuracy
    static const FloatType eps = std::is_same_v<FloatType, float> ? 1e-4 : 1.e-9;
    FloatType totalarea = 0.0;
    for (size_t i = 0; i < smesh.nelements; ++i) {
        REQUIRE(smesh.area(i) > 0);
        REQUIRE(smesh.area(i) == smesh.compute_area(i));
        totalarea += smesh.area(i);
    }
    REQUIRE(std::abs(1 - totalarea / (4.0_ft * pi * EarthRadius * EarthRadius)) < eps);

    // Rotating the pole changes the geometry and must update the cache
    smesh.RotatePoleToGreenland();
    for (size_t i = 0; i < smesh.nelements; ++i) {
        REQUIRE(smesh.area(i) > 0);
        REQUIRE(smesh.area(i) == smesh.compute_area(i));
    }

    // Rotating back must approximately restore the original areas
    // (the rotation round-trip is not bit-exact in floating point)
    smesh.RotatePoleFromGreenland();
    ParametricMesh smesh2(SPHERICAL);
    smesh2.nx = nx;
    smesh2.ny = ny;
    smesh2.nelements = nx * ny;
    smesh2.nnodes = (nx + 1) * (ny + 1);
    smesh2.vertices.resize(smesh2.nnodes, 2);
    for (size_t iy = 0; iy <= ny; ++iy)
        for (size_t ix = 0; ix <= nx; ++ix) {
            smesh2.vertices((nx + 1) * iy + ix, 0) = -1.0_ft * pi + 2.0_ft * pi * ix / nx;
            smesh2.vertices((nx + 1) * iy + ix, 1) = -0.5_ft * pi + 1.0_ft * pi * iy / ny;
        }
    smesh2.landmask.resize(nx * ny, true);
    smesh2.precalculate_area();
    for (size_t i = 0; i < smesh.nelements; ++i)
        REQUIRE(std::abs(1 - smesh.area(i) / smesh2.area(i)) < eps);
}
}
