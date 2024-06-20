/*!
 * @file    XiosGrid_test.cpp
 * @author  Joe Wallwork <jw2423@cam.ac.uk
 * @date    20 June 2024
 * @brief   Tests for XIOS axes
 * @details
 * This test is designed to test axis functionality of the C++ interface
 * for XIOS.
 *
 */
// clang-format off
#include <cstdio>
#include <doctest/extensions/doctest_mpi.h>
#include <iostream>
#include "include/Configurator.hpp"
#include "include/Xios.hpp"
// clang-format on

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
    Nextsim::Configurator::clearStreams();
    std::stringstream config;
    config << "[xios]" << std::endl << "enable = true" << std::endl;
    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Nextsim::Configurator::addStream(std::move(pcstream));

    // Initialize an Xios instance called xios_handler
    Nextsim::Xios xios_handler;
    REQUIRE(xios_handler.isInitialized());

    // Extract MPI size and rank
    const size_t size = xios_handler.getClientMPISize();
    REQUIRE(size == 2);
    const size_t rank = xios_handler.getClientMPIRank();

    // Calendar setup
    xios_handler.setCalendarOrigin(Nextsim::TimePoint("2020-01-23T00:08:15Z"));
    xios_handler.setCalendarStart(Nextsim::TimePoint("2023-03-17T17:11:00Z"));
    xios_handler.setCalendarTimestep(Nextsim::Duration("P0-0T01:30:00"));

    // --- Tests for axis API
    xios_handler.createAxis("axis_A");
    const size_t axis_size = 30;
    xios_handler.setAxisSize("axis_A", axis_size);
    std::vector<double> axisValues(axis_size);
    for (size_t i = 0; i < axis_size; i++) {
        axisValues[i] = i;
    }
    xios_handler.setAxisValues("axis_A", axisValues);

    // --- Tests for domain API
    xios_handler.createDomain("domain_A");
    xios_handler.setDomainType("domain_A", "rectilinear");
    const size_t ni_glo = 60;
    xios_handler.setDomainGlobalLongitudeSize("domain_A", ni_glo);
    const size_t nj_glo = 20;
    xios_handler.setDomainGlobalLatitudeSize("domain_A", nj_glo);
    const size_t ni = ni_glo / size;
    xios_handler.setDomainLongitudeSize("domain_A", ni);
    const size_t nj = nj_glo;
    xios_handler.setDomainLatitudeSize("domain_A", nj);
    xios_handler.setDomainLongitudeStart("domain_A", ni * rank);
    xios_handler.setDomainLatitudeStart("domain_A", 0);
    std::vector<double> vecLon(ni);
    for (size_t i = 0; i < ni; i++) {
        vecLon[i] = -180 + (rank * ni * i) * 360 / ni_glo;
    }
    xios_handler.setDomainLongitudeValues("domain_A", vecLon);
    std::vector<double> vecLat(nj);
    for (size_t j = 0; j < nj; j++) {
        vecLat[j] = -90 + j * 180 / nj_glo;
    }
    xios_handler.setDomainLatitudeValues("domain_A", vecLat);

    // --- Tests for grid API
    const std::string gridId = { "grid_2D" };
    xios_handler.createGrid(gridId);
    // Grid name
    REQUIRE_FALSE(xios_handler.isDefinedGridName(gridId));
    const std::string gridName = { "test_grid" };
    xios_handler.setGridName(gridId, gridName);
    REQUIRE(xios_handler.isDefinedGridName(gridId));
    REQUIRE(xios_handler.getGridName(gridId) == gridName);
    xios_handler.gridAddDomain("grid_2D", "domain_A");
    xios_handler.gridAddAxis("grid_2D", "axis_A");
    // TODO: Test adding domain and axis

    xios_handler.close_context_definition();
    xios_handler.context_finalize();
}
