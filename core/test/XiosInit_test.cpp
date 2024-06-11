/*!
 * @file    XiosInit_test.cpp
 * @author  Joe Wallwork <jw2423@cam.ac.uk
 * @date    7 June 2024
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

    // Initialize an Xios instance called xios_handler
    Nextsim::Xios xios_handler;
    REQUIRE(xios_handler.isInitialized());

    // Extract MPI size and rank
    size_t size = xios_handler.getClientMPISize();
    REQUIRE(size == 2);
    size_t rank = xios_handler.getClientMPIRank();

    // --- Tests for calendar API
    // Calendar type
    std::string calendarType = xios_handler.getCalendarType();
    REQUIRE(calendarType == "Gregorian");
    // Calendar origin
    cxios_date origin;
    origin.year = 2020;
    origin.month = 1;
    origin.day = 23;
    origin.hour = 0;
    origin.minute = 8;
    origin.second = 15;
    std::string datetime = xios_handler.convertXiosDatetimeToString(origin);
    REQUIRE(datetime == "2020-01-23T00:08:15Z");
    xios_handler.setCalendarOrigin(origin);
    datetime = xios_handler.convertXiosDatetimeToString(xios_handler.getCalendarOrigin());
    REQUIRE(datetime == "2020-01-23T00:08:15Z");
    // Calendar start
    cxios_date start;
    start.year = 2023;
    start.month = 3;
    start.day = 17;
    start.hour = 17;
    start.minute = 11;
    start.second = 0;
    datetime = xios_handler.convertXiosDatetimeToString(start);
    REQUIRE(datetime == "2023-03-17T17:11:00Z");
    xios_handler.setCalendarStart(start);
    datetime = xios_handler.convertXiosDatetimeToString(xios_handler.getCalendarStart());
    REQUIRE(datetime == "2023-03-17T17:11:00Z");
    // Timestep
    cxios_duration duration;
    duration.year = 0.0;
    duration.month = 0.0;
    duration.day = 0.0;
    duration.hour = 1.5;
    duration.minute = 0.0;
    duration.second = 0.0;
    duration.timestep = 0.0;
    xios_handler.setCalendarTimestep(duration);
    duration = xios_handler.getCalendarTimestep();
    REQUIRE(duration.year == doctest::Approx(0.0));
    REQUIRE(duration.month == doctest::Approx(0.0));
    REQUIRE(duration.day == doctest::Approx(0.0));
    REQUIRE(duration.hour == doctest::Approx(1.5));
    REQUIRE(duration.minute == doctest::Approx(0.0));
    REQUIRE(duration.second == doctest::Approx(0.0));
    REQUIRE(duration.timestep == doctest::Approx(0.0));

    // --- Tests for axis API
    std::string axisId = { "axis_A" };
    xios_handler.createAxis(axisId);
    // Axis size
    size_t axis_size = 30;
    REQUIRE_FALSE(xios_handler.isDefinedAxisSize(axisId));
    xios_handler.setAxisSize(axisId, axis_size);
    REQUIRE(xios_handler.isDefinedAxisSize(axisId));
    REQUIRE(xios_handler.getAxisSize(axisId) == axis_size);
    // Axis values
    std::vector<double> axisValues;
    for (size_t i = 0; i < axis_size; i++) {
        axisValues.push_back(i);
    }
    REQUIRE_FALSE(xios_handler.areDefinedAxisValues(axisId));
    xios_handler.setAxisValues(axisId, axisValues);
    REQUIRE(xios_handler.areDefinedAxisValues(axisId));
    std::vector<double> axis_A = xios_handler.getAxisValues(axisId);
    for (size_t i = 0; i < axis_size; i++) {
        REQUIRE(axis_A[i] == doctest::Approx(axisValues[i]));
    }

    // --- Tests for domain API
    std::string domainId = { "domain_A" };
    xios_handler.createDomain(domainId);
    // Domain type
    REQUIRE_FALSE(xios_handler.isDefinedDomainType(domainId));
    std::string domainType = { "rectilinear" };
    xios_handler.setDomainType(domainId, domainType);
    REQUIRE(xios_handler.isDefinedDomainType(domainId));
    REQUIRE(xios_handler.getDomainType(domainId) == domainType);
    // Global longitude size
    REQUIRE_FALSE(xios_handler.isDefinedDomainGlobalLongitudeSize(domainId));
    size_t ni_glo = 60;
    xios_handler.setDomainGlobalLongitudeSize(domainId, ni_glo);
    REQUIRE(xios_handler.isDefinedDomainGlobalLongitudeSize(domainId));
    REQUIRE(xios_handler.getDomainGlobalLongitudeSize(domainId) == ni_glo);
    // Global latitude size
    REQUIRE_FALSE(xios_handler.isDefinedDomainGlobalLatitudeSize(domainId));
    size_t nj_glo = 20;
    xios_handler.setDomainGlobalLatitudeSize(domainId, nj_glo);
    REQUIRE(xios_handler.isDefinedDomainGlobalLatitudeSize(domainId));
    REQUIRE(xios_handler.getDomainGlobalLatitudeSize(domainId) == nj_glo);
    // Local longitude size
    REQUIRE_FALSE(xios_handler.isDefinedDomainLongitudeSize(domainId));
    size_t ni = ni_glo / size;
    xios_handler.setDomainLongitudeSize(domainId, ni);
    REQUIRE_FALSE(xios_handler.isDefinedDomainLatitudeSize(domainId));
    REQUIRE(xios_handler.getDomainLongitudeSize(domainId) == ni);
    // Local latitude size
    REQUIRE_FALSE(xios_handler.isDefinedDomainLatitudeSize(domainId));
    size_t nj = nj_glo;
    xios_handler.setDomainLatitudeSize(domainId, nj);
    REQUIRE(xios_handler.isDefinedDomainLatitudeSize(domainId));
    REQUIRE(xios_handler.getDomainLatitudeSize(domainId) == nj);
    // Local longitude start
    REQUIRE_FALSE(xios_handler.isDefinedDomainLongitudeStart(domainId));
    size_t startLon = ni * rank;
    xios_handler.setDomainLongitudeStart(domainId, startLon);
    REQUIRE(xios_handler.isDefinedDomainLongitudeStart(domainId));
    REQUIRE(xios_handler.getDomainLongitudeStart(domainId) == startLon);
    // Local latitude start
    REQUIRE_FALSE(xios_handler.isDefinedDomainLatitudeStart(domainId));
    size_t startLat = 0;
    xios_handler.setDomainLatitudeStart(domainId, startLat);
    REQUIRE(xios_handler.isDefinedDomainLatitudeStart(domainId));
    REQUIRE(xios_handler.getDomainLatitudeStart(domainId) == startLat);
    // Local longitude values
    REQUIRE_FALSE(xios_handler.areDefinedDomainLongitudeValues(domainId));
    std::vector<double> vecLon {};
    for (size_t i = 0; i < ni; i++) {
        vecLon.push_back(-180 + (rank * ni * i) * 360 / ni_glo);
    }
    xios_handler.setDomainLongitudeValues(domainId, vecLon);
    REQUIRE(xios_handler.areDefinedDomainLongitudeValues(domainId));
    std::vector<double> vecLonOut = xios_handler.getDomainLongitudeValues(domainId);
    for (size_t i = 0; i < ni; i++) {
        REQUIRE(vecLonOut[i] == doctest::Approx(vecLon[i]));
    }
    // Local latitude values
    REQUIRE_FALSE(xios_handler.areDefinedDomainLatitudeValues(domainId));
    std::vector<double> vecLat {};
    for (size_t j = 0; j < nj; j++) {
        vecLat.push_back(-90 + j * 180 / nj_glo);
    }
    xios_handler.setDomainLatitudeValues(domainId, vecLat);
    REQUIRE(xios_handler.areDefinedDomainLatitudeValues(domainId));
    std::vector<double> vecLatOut = xios_handler.getDomainLatitudeValues(domainId);
    for (size_t j = 0; j < nj; j++) {
        REQUIRE(vecLatOut[j] == doctest::Approx(vecLat[j]));
    }

    // --- Tests for grid API
    std::string gridId = { "grid_2D" };
    xios_handler.createGrid(gridId);
    // Grid name
    REQUIRE_FALSE(xios_handler.isDefinedGridName(gridId));
    std::string gridName = { "test_grid" };
    xios_handler.setGridName(gridId, gridName);
    REQUIRE(xios_handler.isDefinedGridName(gridId));
    REQUIRE(xios_handler.getGridName(gridId) == gridName);

    // --- Tests for field API
    std::string fieldId = { "field_A" };
    xios_handler.createField(fieldId);
    // Field name
    std::string fieldName = { "test_field" };
    REQUIRE_FALSE(xios_handler.isDefinedFieldName(fieldId));
    xios_handler.setFieldName(fieldId, fieldName);
    REQUIRE(xios_handler.getFieldName(fieldId) == fieldName);
    REQUIRE(xios_handler.isDefinedFieldName(fieldId));
    // Operation
    std::string operation = { "instant" };
    REQUIRE_FALSE(xios_handler.isDefinedFieldOperation(fieldId));
    xios_handler.setFieldOperation(fieldId, operation);
    REQUIRE(xios_handler.isDefinedFieldOperation(fieldId));
    REQUIRE(xios_handler.getFieldOperation(fieldId) == operation);
    // Grid reference
    std::string gridRef = { "grid_2D" };
    REQUIRE_FALSE(xios_handler.isDefinedFieldGridRef(fieldId));
    xios_handler.setFieldGridRef(fieldId, gridRef);
    REQUIRE(xios_handler.isDefinedFieldGridRef(fieldId));
    REQUIRE(xios_handler.getFieldGridRef(fieldId) == gridRef);
    xios_handler.gridAddDomain(gridId, domainId);
    xios_handler.gridAddAxis(gridId, axisId);

    // --- Tests for file API
    std::string fileId { "output" };
    REQUIRE_FALSE(xios_handler.validFileId(fileId));
    xios_handler.createFile(fileId);
    REQUIRE(xios_handler.validFileId(fileId));
    // File name
    std::string fileName { "diagnostic" };
    xios_handler.setFileName(fileId, fileName);
    REQUIRE(xios_handler.isDefinedFileName(fileId));
    REQUIRE(xios_handler.getFileName(fileId) == fileName);
    // File type
    std::string fileType { "one_file" };
    REQUIRE_FALSE(xios_handler.isDefinedFileType(fileId));
    xios_handler.setFileType(fileId, fileType);
    REQUIRE(xios_handler.isDefinedFileType(fileId));
    REQUIRE(xios_handler.getFileType(fileId) == fileType);
    // Output frequency
    REQUIRE_FALSE(xios_handler.isDefinedFileOutputFreq(fileId));
    std::string freq { "1ts" };
    xios_handler.setFileOutputFreq(fileId, freq);
    REQUIRE(xios_handler.isDefinedFileOutputFreq(fileId));
    REQUIRE(xios_handler.getFileOutputFreq(fileId) == freq);
    xios_handler.fileAddField(fileId, fieldId);

    xios_handler.close_context_definition();

    // --- Tests for getCurrentDate method
    REQUIRE(xios_handler.getCurrentDate() == "2023-03-17T17:11:00Z");
    REQUIRE(xios_handler.getCurrentDate(false) == "2023-03-17 17:11:00");

    // create some fake data to test writing methods
    double* field_A = new double[ni * nj * axis_size];
    for (size_t idx = 0; idx < ni * nj * axis_size; idx++) {
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
        xios_handler.write(fieldId, field_A, (int)ni, (int)nj, (int)axis_size);
        // verify timestep
        step = xios_handler.getCalendarStep();
        REQUIRE(step == ts);
    }

    xios_handler.context_finalize();

    // clean up
    delete[] field_A;
}
