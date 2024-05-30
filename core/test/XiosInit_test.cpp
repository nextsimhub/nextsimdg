/*!
 * @file    XiosInit_test.cpp
 * @author  Tom Meltzer <tdm39@cam.ac.uk>
 * @date    Fri 23 Feb 13:43:16 GMT 2024
 * @brief   Tests for XIOS C++ interface
 * @details
 * This test is designed to test all core functionality of the C++ interface
 * for XIOS.
 *
 * Due to the fact that XIOS relies on MPI, it is not convenient to split the
 * tests out into individual test routines. Each one would require a lot of
 * boilerplate, so they have all been group in this test.
 *
 */
#include "include/Configurator.hpp"
#include "include/Xios.hpp"
#include <cstdio>
#include <doctest/extensions/doctest_mpi.h>
#include <iostream>

/*!
 * TestXiosInitialization
 *
 * This test checks all core functionality of the C++ interface for XIOS. It
 * needs to be run with 2 ranks i.e.,
 *
 * `mpirun -n 2 ./testXiosInit_MPI2`
 *
 */
MPI_TEST_CASE("TestXiosInitialization", 2)
{

    // Enable xios in the 'config'
    Nextsim::Configurator::clearStreams();
    std::stringstream config;
    config << "[xios]" << std::endl << "enable = true" << std::endl;
    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Nextsim::Configurator::addStream(std::move(pcstream));

    // initialized instance of xios_handler
    Nextsim::Xios xios_handler;
    REQUIRE(xios_handler.isInitialized());

    // read calendar start date and verify datetime string
    cxios_date start = xios_handler.getCalendarStart();
    std::string datetime = xios_handler.convertXiosDatetimeToString(start);
    REQUIRE(datetime == "2023-03-17T17:11:00Z");

    // read calendar origin date and verify datetime string
    cxios_date origin = xios_handler.getCalendarOrigin();
    datetime = xios_handler.convertXiosDatetimeToString(origin);
    REQUIRE(datetime == "2020-01-23T00:08:15Z");

    // check all elements of cxios_duration struct
    cxios_duration duration;
    duration = xios_handler.getCalendarTimestep();
    REQUIRE(duration.year == doctest::Approx(0.0));
    REQUIRE(duration.month == doctest::Approx(0.0));
    REQUIRE(duration.day == doctest::Approx(0.0));
    REQUIRE(duration.hour == doctest::Approx(1.5));
    REQUIRE(duration.minute == doctest::Approx(0.0));
    REQUIRE(duration.second == doctest::Approx(0.0));
    REQUIRE(duration.timestep == doctest::Approx(0.0));

    // get Calendar start date and modify it
    start = xios_handler.getCalendarStart();
    start.minute = 37;
    // set new start date
    xios_handler.setCalendarStart(start);
    // get Calendar modified start date
    start = xios_handler.getCalendarStart();
    // convert cxios_date to string for comparison
    datetime = xios_handler.convertXiosDatetimeToString(start);
    REQUIRE(datetime == "2023-03-17T17:37:00Z");

    // same steps as calendar start date but for calendar Origin
    origin = xios_handler.getCalendarOrigin();
    origin.second = 1;
    xios_handler.setCalendarOrigin(origin);
    origin = xios_handler.getCalendarOrigin();
    datetime = xios_handler.convertXiosDatetimeToString(origin);
    REQUIRE(datetime == "2020-01-23T00:08:01Z");

    // get Calendar timestep and modify it
    duration = xios_handler.getCalendarTimestep();
    duration.year = 0.5;
    // set Calendar timestep
    xios_handler.setCalendarTimestep(duration);
    // verify timestep has been successfully modified
    duration = xios_handler.getCalendarTimestep();
    REQUIRE(duration.year == doctest::Approx(0.5));
    REQUIRE(duration.month == doctest::Approx(0.0));
    REQUIRE(duration.day == doctest::Approx(0.0));
    REQUIRE(duration.hour == doctest::Approx(1.5));
    REQUIRE(duration.minute == doctest::Approx(0.0));
    REQUIRE(duration.second == doctest::Approx(0.0));
    REQUIRE(duration.timestep == doctest::Approx(0.0));

    xios_handler.context_finalize();

    // check the getCurrentDate method
    xios::CDate current;
    current = xios_handler.getCurrentDate();
    std::string current_date = {};
    current_date = current.toString().c_str();
    REQUIRE(current_date == "2023-03-17 17:37:00");

    // create some fake data to test writing methods
    int ni = 30;
    int nj = 30;
    double* field_A = new double[ni * nj];
    for (int idx = 0; idx < ni * nj; idx++) {
        field_A[idx] = 1.0 * idx;
    }

    // create field
    std::string fieldId = { "field_A" };

    // verify calendar step is starting from zero
    int step = xios_handler.getCalendarStep();
    REQUIRE(step == 0);

    // simulate 4 iterations (timesteps)
    for (int ts = 1; ts <= 4; ts++) {
        // update the current timestep
        xios_handler.updateCalendar(ts);
        // send data to XIOS to be written to disk
        xios_handler.write(fieldId, field_A, ni, nj);
        // verify timestep
        step = xios_handler.getCalendarStep();
        REQUIRE(step == ts);
    }

    // clean up
    delete[] field_A;
}
