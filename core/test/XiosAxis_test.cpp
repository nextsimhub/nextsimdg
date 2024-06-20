/*!
 * @file    XiosAxis_test.cpp
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
    Nextsim::Configurator::clearStreams();
    std::stringstream config;
    config << "[xios]" << std::endl << "enable = true" << std::endl;
    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Nextsim::Configurator::addStream(std::move(pcstream));

    // Initialize an Xios instance called xios_handler
    Nextsim::Xios xios_handler;
    REQUIRE(xios_handler.isInitialized());
    REQUIRE(xios_handler.getClientMPISize() == 2);

    // Calendar setup
    xios_handler.setCalendarOrigin(Nextsim::TimePoint("2020-01-23T00:08:15Z"));
    xios_handler.setCalendarStart(Nextsim::TimePoint("2023-03-17T17:11:00Z"));
    xios_handler.setCalendarTimestep(Nextsim::Duration("P0-0T01:30:00"));

    // --- Tests for axis API
    const std::string axisId = { "axis_A" };
    xios_handler.createAxis(axisId);
    // Axis size
    const size_t axis_size = 30;
    REQUIRE_FALSE(xios_handler.isDefinedAxisSize(axisId));
    xios_handler.setAxisSize(axisId, axis_size);
    REQUIRE(xios_handler.isDefinedAxisSize(axisId));
    REQUIRE(xios_handler.getAxisSize(axisId) == axis_size);
    // Axis values
    std::vector<double> axisValues(axis_size);
    for (size_t i = 0; i < axis_size; i++) {
        axisValues[i] = i;
    }
    REQUIRE_FALSE(xios_handler.areDefinedAxisValues(axisId));
    xios_handler.setAxisValues(axisId, axisValues);
    REQUIRE(xios_handler.areDefinedAxisValues(axisId));
    std::vector<double> axis_A = xios_handler.getAxisValues(axisId);
    for (size_t i = 0; i < axis_size; i++) {
        REQUIRE(axis_A[i] == doctest::Approx(axisValues[i]));
    }

    xios_handler.close_context_definition();
    xios_handler.context_finalize();
}
