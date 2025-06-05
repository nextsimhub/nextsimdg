/*!
 * @file ParametricMesh_test.cpp
 *
 * @brief Test the ParametricMesh class, especially processing from ModelArray files.
 *
 * @date 09 Jan 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "include/ParametricMesh.hpp"

#include "../test/FakeSmeshData.hpp"
#include "include/gridNames.hpp"

#include <filesystem>

namespace Nextsim {
#define TO_STR(s) TO_STRI(s)
#define TO_STRI(s) #s
#ifndef TEST_FILE_SOURCE
#define TEST_FILE_SOURCE .
#endif

static const size_t nx = 154;
static const size_t ny = 121;
static const double dx = 25000;
static const double dy = 25000;

TEST_SUITE_BEGIN("ParametricMeshArea");

TEST_CASE("Test area cartesian")
{
    ParametricMesh smesh(CARTESIAN);

    // define mesh
    smesh.nx = nx;
    smesh.ny = ny;
    smesh.nelements = nx * ny;
    smesh.nnodes = (nx + 1) * (ny + 1);
    smesh.vertices.resize(smesh.nnodes, 2);
    for (int iy = 0; iy <= ny; ++iy)
        for (int ix = 0; ix <= nx; ++ix) {
            smesh.vertices((nx + 1) * iy + ix, 0) = dx * ix / nx;
            smesh.vertices((nx + 1) * iy + ix, 1) = dy * iy / ny;
        }
    smesh.landmask.resize(nx * ny, true);

    // exact area
    double exact = dx * dy;

    // check mesh area
    double totalarea = 0.0;
    for (int i = 0; i < nx * ny; ++i) {
        REQUIRE(smesh.area(i) > 0);
        totalarea += smesh.area(i); // area must be positive
    }

    // total area must match domain dimension up to machine accuracy
    REQUIRE(fabs(1.0 - totalarea / exact) < 1.e-12);

    //////// next, distort the inner vertices of the mesh
    double h = std::min(dx / nx, dy / ny); // min. mesh size in x/y-direction
    for (int iy = 1; iy < ny; ++iy)
        for (int ix = 1; ix < nx; ++ix) {
            smesh.vertices((nx + 1) * iy + ix, 0) += 0.2 * h * sin(ix + iy); // distort by 20%
            smesh.vertices((nx + 1) * iy + ix, 1) += 0.2 * h * cos(ix + iy);
        }

    // check mesh area
    totalarea = 0.0;
    for (int i = 0; i < nx * ny; ++i) {
        REQUIRE(smesh.area(i) > 0);
        totalarea += smesh.area(i); // area must be positive
    }
    // total area must match domain dimension up to machine accuracy
    REQUIRE(fabs(1.0 - totalarea / exact) < 1.e-12);
}

TEST_CASE("Test area spherical")
{
    ParametricMesh smesh(SPHERICAL);

    // define mesh
    smesh.nx = nx;
    smesh.ny = ny;
    smesh.nelements = nx * ny;
    smesh.nnodes = (nx + 1) * (ny + 1);
    smesh.vertices.resize(smesh.nnodes, 2);
    for (int iy = 0; iy <= ny; ++iy)
        for (int ix = 0; ix <= nx; ++ix) {
            smesh.vertices((nx + 1) * iy + ix, 0) = -1.0 * M_PI + 2.0 * M_PI * ix / nx;
            smesh.vertices((nx + 1) * iy + ix, 1) = -0.5 * M_PI + 1.0 * M_PI * iy / ny;
        }
    smesh.landmask.resize(nx * ny, true);

    // exact radius (surface of earth)
    double exact = 4.0 * M_PI * EarthRadius * EarthRadius;

    // check mesh area
    double totalarea = 0.0;
    for (int i = 0; i < nx * ny; ++i) {
        REQUIRE(smesh.area(i) > 0);
        totalarea += smesh.area(i); // area must be positive
    }

    // total area must match domain dimension up to machine accuracy
    REQUIRE(fabs(1.0 - totalarea / exact) < 1.e-9);

    //////// next, distort the inner vertices of the mesh
    double h = std::min(2.0 * M_PI / nx, M_PI / ny); // min. mesh size in x/y-direction
    for (int iy = 1; iy < ny; ++iy)
        for (int ix = 1; ix < nx; ++ix) {
            smesh.vertices((nx + 1) * iy + ix, 0) += 0.2 * h * sin(ix + iy); // distort by 20%
            smesh.vertices((nx + 1) * iy + ix, 1) += 0.2 * h * cos(ix + iy);
        }

    // check mesh area
    totalarea = 0.0;
    for (int i = 0; i < nx * ny; ++i) {
        REQUIRE(smesh.area(i) > 0);
        totalarea += smesh.area(i); // area must be positive
    }
    // total area must match domain dimension up to machine accuracy
    REQUIRE(fabs(1.0 - totalarea / exact) < 1.e-9);
}
}
