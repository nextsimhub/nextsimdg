/*!
 * @file    XiosAxis_test.cpp
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @date    12 August 2024
 * @brief   Tests for XIOS axes
 * @details
 * This test is designed to test axis functionality of the C++ interface
 * for XIOS.
 *
 */
#include <doctest/extensions/doctest_mpi.h>
#undef INFO

#include "StructureModule/include/ParametricGrid.hpp"
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
    REQUIRE_THROWS_WITH(xios_handler.getAxisSize(axisId), "Xios: Undefined axis 'axis_A'");
    REQUIRE_THROWS_WITH(xios_handler.getAxisValues(axisId), "Xios: Undefined axis 'axis_A'");
    xios_handler.createAxis(axisId);
    REQUIRE_THROWS_WITH(xios_handler.createAxis(axisId), "Xios: Axis 'axis_A' already exists");
    // Axis size
    REQUIRE_THROWS_WITH(xios_handler.getAxisSize(axisId), "Xios: Undefined size for axis 'axis_A'");
    const size_t axisSize { 2 };
    xios_handler.setAxisSize(axisId, axisSize);
    REQUIRE(xios_handler.getAxisSize(axisId) == axisSize);
    // Axis values
    REQUIRE_THROWS_WITH(
        xios_handler.getAxisValues(axisId), "Xios: Undefined values for axis 'axis_A'");
    REQUIRE_THROWS_WITH(xios_handler.setAxisValues(axisId, { 0.0, 1.0, 2.0 }),
        "Xios: Size incompatible with values for axis 'axis_A'");
    std::vector<double> axisValues { 0.0, 1.0 };
    xios_handler.setAxisValues(axisId, axisValues);
    std::vector<double> axis_A = xios_handler.getAxisValues(axisId);
    REQUIRE(axis_A[0] == doctest::Approx(0.0));
    REQUIRE(axis_A[1] == doctest::Approx(1.0));

    xios_handler.close_context_definition();
    xios_handler.context_finalize();
}
}
