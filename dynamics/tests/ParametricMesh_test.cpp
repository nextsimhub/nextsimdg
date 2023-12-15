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

#include "FakeSmeshData.hpp"

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

    // Dirichlet conditions
    // Element 5 comes first, and is closed on edges 0 & 3
    REQUIRE(smesh.dirichlet[0][0] == 5);
    REQUIRE(smesh.dirichlet[3][0] == 5);
    // Element 7 is the first with a closed edge 1
    REQUIRE(smesh.dirichlet[1][0] == 7);
    // Element 48 is the first with a closed edge 2
    REQUIRE(smesh.dirichlet[2][0] == 48);
    // sizes of the Dirichlet arrays
    REQUIRE(smesh.dirichlet[0].size() == 369);
    REQUIRE(smesh.dirichlet[1].size() == 384);
    REQUIRE(smesh.dirichlet[0].size() == smesh.dirichlet[2].size());
    REQUIRE(smesh.dirichlet[1].size() == smesh.dirichlet[3].size());

    // No periodic boundary conditions
    REQUIRE(smesh.periodic.size() == 0);
}

TEST_CASE("Compare readmesh and landmask reading")
{

    ParametricMesh fromFile(CoordinateSystem);
    ParametricMesh fromArrays(CoordinateSystem);

    fromFile.readmesh(smeshFile);


}
}
