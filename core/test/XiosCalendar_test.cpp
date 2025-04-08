/*!
 * @file    XiosCalendar_test.cpp
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @date    08 Apr 2025
 * @brief   Tests for XIOS calandars
 * @details
 * This test is designed to test calendar functionality of the C++ interface
 * for XIOS.
 *
 */
#include <doctest/extensions/doctest_mpi.h>
#undef INFO

#include "StructureModule/include/ParametricGrid.hpp"
#include "include/Xios.hpp"

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
    enableXios();

    // Get the Xios singleton instance and check it's initialized
    Xios& xiosHandler = Xios::getInstance();
    REQUIRE(xiosHandler.isInitialized());
    REQUIRE(xiosHandler.getClientMPISize() == 2);

    // --- Tests for calendar API
    // Calendar type
    REQUIRE(xiosHandler.getCalendarType() == "Gregorian");
    // Calendar origin
    REQUIRE(xiosHandler.getCalendarOrigin().format() == "1970-01-01T00:00:00Z"); // Default
    TimePoint origin("2020-01-23T00:08:15Z");
    xiosHandler.setCalendarOrigin(origin);
    REQUIRE(origin == xiosHandler.getCalendarOrigin());
    REQUIRE(origin.format() == "2020-01-23T00:08:15Z");
    // Calendar start
    REQUIRE(xiosHandler.getCalendarStart().format() == "1970-01-01T00:00:00Z"); // Default
    TimePoint start("2023-03-17T17:11:00Z");
    xiosHandler.setCalendarStart(start);
    REQUIRE(start == xiosHandler.getCalendarStart());
    REQUIRE(start.format() == "2023-03-17T17:11:00Z");
    // Timestep
    REQUIRE(xiosHandler.getCalendarTimestep().seconds() == 3600.0); // Default
    Duration timestep("P0-0T01:30:00");
    REQUIRE(timestep.seconds() == 5400.0);
    xiosHandler.setCalendarTimestep(timestep);
    REQUIRE(xiosHandler.getCalendarTimestep().seconds() == 5400.0);

    xiosHandler.close_context_definition();

    // --- Tests for getCurrentDate method
    REQUIRE(xiosHandler.getCalendarStep() == 0);
    REQUIRE(xiosHandler.getCurrentDate() == "2023-03-17T17:11:00Z");
    REQUIRE(xiosHandler.getCurrentDate(false) == "2023-03-17 17:11:00");

    // -- Tests that the timestep is set up correctly
    xiosHandler.setCalendarStep(1);
    REQUIRE(xiosHandler.getCurrentDate() == "2023-03-17T18:41:00Z");

    xiosHandler.context_finalize();
}

}
