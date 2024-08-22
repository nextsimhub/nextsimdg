/*!
 * @file    XiosGrid_test.cpp
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @date    5 August 2024
 * @brief   Tests for XIOS axes
 * @details
 * This test is designed to test axis functionality of the C++ interface
 * for XIOS.
 *
 */
#include <doctest/extensions/doctest_mpi.h>
#undef INFO

#include "include/Configurator.hpp"
#include "include/Xios.hpp"

#include <iostream>

namespace Nextsim {

/*!
 * TestXiosGrid
 *
 * This function tests the grid functionality of the C++ interface for XIOS. It
 * needs to be run with 2 ranks i.e.,
 *
 * `mpirun -n 2 ./testXiosGrid_MPI2`
 *
 */
MPI_TEST_CASE("TestXiosGrid", 2)
{

    // Enable XIOS in the 'config'
    Configurator::clearStreams();
    std::stringstream config;
    config << "[xios]" << std::endl << "enable = true" << std::endl;
    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    // Initialize an Xios instance called xios_handler
    Xios xios_handler;
    REQUIRE(xios_handler.isInitialized());
    const size_t size = xios_handler.getClientMPISize();
    REQUIRE(size == 2);
    const size_t rank = xios_handler.getClientMPIRank();

    // Set timestep as a minimum
    xios_handler.setCalendarTimestep(Duration("P0-0T01:30:00"));

    // Create a 4x2 horizontal domain with a partition halving the x-extent
    xios_handler.createDomain("xy_domain");
    xios_handler.setDomainType("xy_domain", "rectilinear");
    xios_handler.setDomainGlobalXSize("xy_domain", 4);
    xios_handler.setDomainGlobalYSize("xy_domain", 2);
    xios_handler.setDomainLocalXStart("xy_domain", 2 * rank);
    xios_handler.setDomainLocalYStart("xy_domain", 0);
    xios_handler.setDomainLocalXValues("xy_domain", { -1.0 + rank, -0.5 + rank });
    xios_handler.setDomainLocalYValues("xy_domain", { -1.0, 1.0 });

    // Create a vertical axis with 2 points
    xios_handler.createAxis("z_axis");
    xios_handler.setAxisValues("z_axis", { 0.0, 1.0 });

    // --- Tests for grid API
    const std::string gridId = { "grid_2D" };
    REQUIRE_THROWS_WITH(xios_handler.getGridName(gridId), "Xios: Undefined grid 'grid_2D'");
    xios_handler.createGrid(gridId);
    REQUIRE_THROWS_WITH(xios_handler.createGrid(gridId), "Xios: Grid 'grid_2D' already exists");
    // Grid name
    const std::string gridName = { "test_grid" };
    REQUIRE_THROWS_WITH(
        xios_handler.getGridName(gridId), "Xios: Undefined name for grid 'grid_2D'");
    xios_handler.setGridName(gridId, gridName);
    REQUIRE(xios_handler.getGridName(gridId) == gridName);
    // Add axis
    xios_handler.gridAddAxis("grid_2D", "z_axis");
    std::vector<std::string> axisIds = xios_handler.gridGetAxisIds(gridId);
    REQUIRE(axisIds.size() == 1);
    REQUIRE(axisIds[0] == "z_axis");
    // Add domain
    xios_handler.gridAddDomain("grid_2D", "xy_domain");
    std::vector<std::string> domainIds = xios_handler.gridGetDomainIds(gridId);
    REQUIRE(domainIds.size() == 1);
    REQUIRE(domainIds[0] == "xy_domain");

    xios_handler.close_context_definition();
    xios_handler.context_finalize();
}
}
