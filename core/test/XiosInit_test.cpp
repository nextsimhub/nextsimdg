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
// clang-format off
#include <cstdio>
#include <doctest/extensions/doctest_mpi.h>
#include <iostream>
#include "include/Configurator.hpp"
#include "include/Xios.hpp"
// clang-format on

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

    xios_handler.close_context_definition();

    // check the getCurrentDate method with and without ISO formatting
    xios::CDate current;
    std::string current_date = xios_handler.getCurrentDate();
    REQUIRE(current_date == "2023-03-17T17:37:00Z");
    current_date = xios_handler.getCurrentDate(false);
    REQUIRE(current_date == "2023-03-17 17:37:00");

    // check axis getters
    int axis_size = xios_handler.getAxisSize("axis_A");
    REQUIRE(axis_size == 30);
    std::vector<double> axis_A = xios_handler.getAxisValues("axis_A");
    for (int i = 0; i < axis_size; i++) {
        REQUIRE(axis_A[i] == doctest::Approx(i));
    }

    // check global domain getters
    std::string domainId { "domain_A" };
    REQUIRE(xios_handler.getDomainType(domainId) == "rectilinear");
    int ni_glo = xios_handler.getDomainGlobalLongitudeSize(domainId);
    int nj_glo = xios_handler.getDomainGlobalLatitudeSize(domainId);
    REQUIRE(ni_glo == 40);
    REQUIRE(nj_glo == 60);

    // check local domain setters and getters
    int rank = xios_handler.rank;
    int ni = ni_glo / xios_handler.size;
    int nj = nj_glo;
    xios_handler.setDomainLongitudeSize(domainId, ni);
    xios_handler.setDomainLatitudeSize(domainId, nj);
    REQUIRE(xios_handler.getDomainLongitudeSize(domainId) == ni);
    REQUIRE(xios_handler.getDomainLatitudeSize(domainId) == nj);
    int startLon = ni * rank;
    int startLat = 0;
    xios_handler.setDomainLongitudeStart(domainId, startLon);
    xios_handler.setDomainLatitudeStart(domainId, startLat);
    REQUIRE(xios_handler.getDomainLongitudeStart(domainId) == startLon);
    REQUIRE(xios_handler.getDomainLatitudeStart(domainId) == startLat);
    std::vector<double> vecLon {};
    std::vector<double> vecLat {};
    for (int i = 0; i < ni; i++) {
        vecLon.push_back(-180 + (rank * ni * i) * 360 / ni_glo);
    }
    for (int j = 0; j < nj; j++) {
        vecLat.push_back(-90 + j * 180 / nj_glo);
    }
    xios_handler.setDomainLongitudeValues(domainId, vecLon);
    xios_handler.setDomainLatitudeValues(domainId, vecLat);
    std::vector<double> vecLonOut = xios_handler.getDomainLongitudeValues(domainId);
    std::vector<double> vecLatOut = xios_handler.getDomainLatitudeValues(domainId);
    for (int i = 0; i < ni; i++) {
        REQUIRE(vecLonOut[i] == doctest::Approx(vecLon[i]));
    }
    for (int j = 0; j < nj; j++) {
        REQUIRE(vecLatOut[j] == doctest::Approx(vecLat[j]));
    }

    // check field getters
    std::string fieldId = { "field_A" };
    REQUIRE(xios_handler.isDefinedFieldName(fieldId));
    REQUIRE(xios_handler.isDefinedFieldOperation(fieldId));
    REQUIRE(xios_handler.isDefinedFieldGridRef(fieldId));
    REQUIRE(xios_handler.getFieldName(fieldId) == "test_field");
    REQUIRE(xios_handler.getFieldOperation(fieldId) == "instant");
    REQUIRE(xios_handler.getFieldGridRef(fieldId) == "grid_2D");

    // check grid getters
    std::string gridId = { "grid_2D" };
    REQUIRE(xios_handler.getGridName(gridId) == "test_grid");

    // check file getters
    REQUIRE_FALSE(xios_handler.validFileId("invalid"));
    std::string fileId { "output" };
    REQUIRE(xios_handler.validFileId(fileId));
    REQUIRE(xios_handler.getFileName(fileId) == "diagnostic");
    REQUIRE(xios_handler.getFileType(fileId) == "one_file");
    REQUIRE(xios_handler.isDefinedFileOutputFreq(fileId));
    REQUIRE(xios_handler.getFileOutputFreq(fileId) == "1ts");

    // create some fake data to test writing methods
    double* field_A = new double[ni * nj * axis_size];
    for (int idx = 0; idx < ni * nj * axis_size; idx++) {
        field_A[idx] = 1.0 * idx;
    }

    // verify calendar step is starting from zero
    int step = xios_handler.getCalendarStep();
    REQUIRE(step == 0);

    // simulate 4 iterations (timesteps)
    for (int ts = 1; ts <= 4; ts++) {
        // update the current timestep
        xios_handler.updateCalendar(ts);
        // send data to XIOS to be written to disk
        // xios_handler.write(fieldId, field_A, ni, nj, axis_size); // FIXME
        // verify timestep
        step = xios_handler.getCalendarStep();
        REQUIRE(step == ts);
    }

    xios_handler.context_finalize();

    // clean up
    delete[] field_A;
}
