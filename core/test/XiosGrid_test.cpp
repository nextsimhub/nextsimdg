/*!
 * @file    XiosGrid_test.cpp
 * @author  Joe Wallwork <jw2423@cam.ac.uk
 * @date    26 July 2024
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

    // Axis setup
    xios_handler.createAxis("axis_A");
    xios_handler.setAxisValues("axis_A", { 0, 1 });

    // Domain setup
    xios_handler.createDomain("domain_A");
    xios_handler.setDomainType("domain_A", "rectilinear");
    const size_t nx_glo = 4;
    xios_handler.setDomainGlobalXSize("domain_A", nx_glo);
    const size_t ny_glo = 2;
    xios_handler.setDomainGlobalYSize("domain_A", ny_glo);
    const size_t nx = nx_glo / size;
    xios_handler.setDomainLocalXSize("domain_A", nx);
    const size_t ny = ny_glo;
    xios_handler.setDomainLocalYSize("domain_A", ny);
    xios_handler.setDomainLocalXStart("domain_A", nx * rank);
    xios_handler.setDomainLocalYStart("domain_A", 0);
    std::vector<double> vx { -1.0 + rank, -0.5 + rank };
    xios_handler.setDomainLocalXValues("domain_A", vx);
    std::vector<double> vy { -1.0, 1.0 };
    xios_handler.setDomainLocalYValues("domain_A", vy);

    // --- Tests for grid API
    const std::string gridId = { "grid_2D" };
    xios_handler.createGrid(gridId);
    // Grid name
    const std::string gridName = { "test_grid" };
    xios_handler.setGridName(gridId, gridName);
    REQUIRE(xios_handler.getGridName(gridId) == gridName);
    // Add axis
    xios_handler.gridAddAxis("grid_2D", "axis_A");
    std::vector<std::string> axisIds = xios_handler.gridGetAxisIds(gridId);
    REQUIRE(axisIds.size() == 1);
    REQUIRE(axisIds[0] == "axis_A");
    // Add domain
    xios_handler.gridAddDomain("grid_2D", "domain_A");
    std::vector<std::string> domainIds = xios_handler.gridGetDomainIds(gridId);
    REQUIRE(domainIds.size() == 1);
    REQUIRE(domainIds[0] == "domain_A");

    xios_handler.close_context_definition();
    xios_handler.context_finalize();
}
}
