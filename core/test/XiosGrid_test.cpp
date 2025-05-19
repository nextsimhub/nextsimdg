/*!
 * @file    XiosGrid_test.cpp
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @author  Adeleke Bankole <ab3191@cam.ac.uk>
 * @date    23 Apr 2025
 * @brief   Tests for XIOS grid
 * @details
 * This test is designed to test grid functionality of the C++ interface
 * for XIOS.
 *
 */
#include <doctest/extensions/doctest_mpi.h>
#undef INFO

#include "include/Finalizer.hpp"
#include "include/Xios.hpp"

namespace Nextsim {

/*!
 * TestXiosGrid
 *
 * This function tests the grid functionality of the C++ interface for XIOS. It
 * needs to be run with 4 ranks i.e.,
 *
 * `mpirun -n 4 ./testXiosGrid_MPI4`
 *
 */
MPI_TEST_CASE("TestXiosGrid", 4)
{
    enableXios();
    // Get the Xios singleton instance and check it's initialized
    Xios& xiosHandler = Xios::getInstance();
    REQUIRE(xiosHandler.isInitialized());
    const size_t size = xiosHandler.getClientMPISize();
    REQUIRE(size == 4);
    const size_t rank = xiosHandler.getClientMPIRank();

    // Create a 4x2 horizontal domain with a partition halving the x-extent
    xiosHandler.createDomain("domain_XY");
    xiosHandler.setDomainType("domain_XY", "rectilinear");
    xiosHandler.setDomainGlobalXSize("domain_XY", 4);
    xiosHandler.setDomainGlobalYSize("domain_XY", 2);
    xiosHandler.setDomainLocalXStart("domain_XY", 2 * rank);
    xiosHandler.setDomainLocalYStart("domain_XY", 0);
    xiosHandler.setDomainLocalXValues("domain_XY", { -1.0 + rank, -0.5 + rank });
    xiosHandler.setDomainLocalYValues("domain_XY", { -1.0, 1.0 });

    // Create a vertical axis with 2 points
    xiosHandler.createAxis("axis_Z");
    xiosHandler.setAxisValues("axis_Z", { 0.0, 1.0 });

    // --- Tests for grid API
    const std::string gridId = { "grid_2D" };
    REQUIRE_THROWS_WITH(xiosHandler.getGridName(gridId), "Xios: Undefined grid 'grid_2D'");
    xiosHandler.createGrid(gridId);
    REQUIRE_THROWS_WITH(xiosHandler.createGrid(gridId), "Xios: Grid 'grid_2D' already exists");
    // Grid name
    const std::string gridName = { "test_grid" };
    REQUIRE_THROWS_WITH(xiosHandler.getGridName(gridId), "Xios: Undefined name for grid 'grid_2D'");
    xiosHandler.setGridName(gridId, gridName);
    REQUIRE(xiosHandler.getGridName(gridId) == gridName);
    // Add axis
    xiosHandler.gridAddAxis("grid_2D", "axis_Z");
    std::vector<std::string> axisIds = xiosHandler.gridGetAxisIds(gridId);
    REQUIRE(axisIds.size() == 1);
    REQUIRE(axisIds[0] == "axis_Z");
    // Add domain
    xiosHandler.gridAddDomain("grid_2D", "domain_XY");
    std::vector<std::string> domainIds = xiosHandler.gridGetDomainIds(gridId);
    REQUIRE(domainIds.size() == 1);
    REQUIRE(domainIds[0] == "domain_XY");

    xiosHandler.close_context_definition();
    xiosHandler.context_finalize();
    Finalizer::finalize();
}
}
