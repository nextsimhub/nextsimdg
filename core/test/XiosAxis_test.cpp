/*!
 * @file    XiosAxis_test.cpp
 * @author  Joe Wallwork <jw2423@cam.ac.uk
 * @date    27 June 2024
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
 * TestXiosAxis
 *
 * This function tests the axis functionality of the C++ interface for XIOS. It
 * needs to be run with 2 ranks i.e.,
 *
 * `mpirun -n 2 ./testXiosAxis_MPI2`
 *
 */
MPI_TEST_CASE("TestXiosAxis", 2)
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
    REQUIRE(xios_handler.getClientMPISize() == 2);

    // Set timestep as a minimum
    xios_handler.setCalendarTimestep(Duration("P0-0T01:00:00"));

    // --- Tests for axis API
    const std::string axisId = { "axis_A" };
    xios_handler.createAxis(axisId);
    // Axis size
    const size_t axis_size { 2 };
    REQUIRE_FALSE(xios_handler.isDefinedAxisSize(axisId));
    xios_handler.setAxisSize(axisId, axis_size);
    REQUIRE(xios_handler.isDefinedAxisSize(axisId));
    REQUIRE(xios_handler.getAxisSize(axisId) == axis_size);
    // Axis values
    std::vector<double> axisValues { 0, 1 };
    REQUIRE_FALSE(xios_handler.areDefinedAxisValues(axisId));
    xios_handler.setAxisValues(axisId, axisValues);
    REQUIRE(xios_handler.areDefinedAxisValues(axisId));
    std::vector<double> axis_A = xios_handler.getAxisValues(axisId);
    REQUIRE(axis_A[0] == doctest::Approx(0));
    REQUIRE(axis_A[1] == doctest::Approx(1));

    xios_handler.close_context_definition();
    xios_handler.context_finalize();
}
}
