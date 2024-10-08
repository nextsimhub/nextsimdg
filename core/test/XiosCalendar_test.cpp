/*!
 * @file    XiosCalendar_test.cpp
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @date    5 August 2024
 * @brief   Tests for XIOS calandars
 * @details
 * This test is designed to test calendar functionality of the C++ interface
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
 * TestXiosCalendar
 *
 * This function tests the calendar functionality of the C++ interface for XIOS. It
 * needs to be run with 2 ranks i.e.,
 *
 * `mpirun -n 2 ./testXiosCalendar_MPI2`
 *
 */
MPI_TEST_CASE("TestXiosInitialization", 2)
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

    // --- Tests for calendar API
    // Calendar type
    REQUIRE(xios_handler.getCalendarType() == "Gregorian");
    // Calendar origin
    TimePoint origin("2020-01-23T00:08:15Z");
    xios_handler.setCalendarOrigin(origin);
    REQUIRE(origin == xios_handler.getCalendarOrigin());
    REQUIRE(origin.format() == "2020-01-23T00:08:15Z");
    // Calendar start
    TimePoint start("2023-03-17T17:11:00Z");
    xios_handler.setCalendarStart(start);
    REQUIRE(start == xios_handler.getCalendarStart());
    REQUIRE(start.format() == "2023-03-17T17:11:00Z");
    // Timestep
    Duration timestep("P0-0T01:30:00");
    REQUIRE(timestep.seconds() == doctest::Approx(5400.0));
    xios_handler.setCalendarTimestep(timestep);
    REQUIRE(xios_handler.getCalendarTimestep().seconds() == doctest::Approx(5400.0));

    xios_handler.close_context_definition();

    // --- Tests for getCurrentDate method
    REQUIRE(xios_handler.getCalendarStep() == 0);
    REQUIRE(xios_handler.getCurrentDate() == "2023-03-17T17:11:00Z");
    REQUIRE(xios_handler.getCurrentDate(false) == "2023-03-17 17:11:00");

    // -- Tests that the timestep is set up correctly
    xios_handler.updateCalendar(1);
    REQUIRE(xios_handler.getCurrentDate() == "2023-03-17T18:41:00Z");

    xios_handler.context_finalize();
}

}
