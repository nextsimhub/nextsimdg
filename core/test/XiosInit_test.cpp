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

    // --- Tests for axis API
    std::string axisId = { "axis_A" };
    // xios_handler.createAxis(axisId); // FIXME
    // Axis size
    int axis_size = 30;
    REQUIRE_FALSE(xios_handler.isDefinedAxisSize(axisId));
    xios_handler.setAxisSize(axisId, axis_size);
    REQUIRE(xios_handler.isDefinedAxisSize(axisId));
    REQUIRE(xios_handler.getAxisSize(axisId) == axis_size);
    // Axis values
    std::vector<double> axisValues;
    for (int i = 0; i < axis_size; i++) {
        axisValues.push_back(i);
    }
    REQUIRE_FALSE(xios_handler.areDefinedAxisValues(axisId));
    xios_handler.setAxisValues(axisId, axisValues);
    REQUIRE(xios_handler.areDefinedAxisValues(axisId));
    std::vector<double> axis_A = xios_handler.getAxisValues(axisId);
    for (int i = 0; i < axis_size; i++) {
        REQUIRE(axis_A[i] == doctest::Approx(axisValues[i]));
    }

    // --- Tests for domain API
    std::string domainId = { "domain_A" };
    // xios_handler.createDomain(domainId); // FIXME
    // Domain type
    REQUIRE_FALSE(xios_handler.isDefinedDomainType(domainId));
    std::string domainType = { "rectilinear" };
    xios_handler.setDomainType(domainId, domainType);
    REQUIRE(xios_handler.isDefinedDomainType(domainId));
    REQUIRE(xios_handler.getDomainType(domainId) == domainType);
    // Global longitude size
    REQUIRE_FALSE(xios_handler.isDefinedDomainGlobalLongitudeSize(domainId));
    int ni_glo = 60;
    xios_handler.setDomainGlobalLongitudeSize(domainId, ni_glo);
    REQUIRE(xios_handler.isDefinedDomainGlobalLongitudeSize(domainId));
    REQUIRE(xios_handler.getDomainGlobalLongitudeSize(domainId) == ni_glo);
    // Global latitude size
    REQUIRE_FALSE(xios_handler.isDefinedDomainGlobalLatitudeSize(domainId));
    int nj_glo = 20;
    xios_handler.setDomainGlobalLatitudeSize(domainId, nj_glo);
    REQUIRE(xios_handler.isDefinedDomainGlobalLatitudeSize(domainId));
    REQUIRE(xios_handler.getDomainGlobalLatitudeSize(domainId) == nj_glo);
    // Local longitude size
    REQUIRE_FALSE(xios_handler.isDefinedDomainLongitudeSize(domainId));
    int ni = ni_glo / xios_handler.size;
    xios_handler.setDomainLongitudeSize(domainId, ni);
    REQUIRE_FALSE(xios_handler.isDefinedDomainLatitudeSize(domainId));
    REQUIRE(xios_handler.getDomainLongitudeSize(domainId) == ni);
    // Local latitude size
    REQUIRE_FALSE(xios_handler.isDefinedDomainLatitudeSize(domainId));
    int nj = nj_glo;
    xios_handler.setDomainLatitudeSize(domainId, nj);
    REQUIRE(xios_handler.isDefinedDomainLatitudeSize(domainId));
    REQUIRE(xios_handler.getDomainLatitudeSize(domainId) == nj);
    // Local longitude start
    REQUIRE_FALSE(xios_handler.isDefinedDomainLongitudeStart(domainId));
    int rank = xios_handler.rank;
    int startLon = ni * rank;
    xios_handler.setDomainLongitudeStart(domainId, startLon);
    REQUIRE(xios_handler.isDefinedDomainLongitudeStart(domainId));
    REQUIRE(xios_handler.getDomainLongitudeStart(domainId) == startLon);
    // Local latitude start
    REQUIRE_FALSE(xios_handler.isDefinedDomainLatitudeStart(domainId));
    int startLat = 0;
    xios_handler.setDomainLatitudeStart(domainId, startLat);
    REQUIRE(xios_handler.isDefinedDomainLatitudeStart(domainId));
    REQUIRE(xios_handler.getDomainLatitudeStart(domainId) == startLat);
    // Local longitude values
    REQUIRE_FALSE(xios_handler.areDefinedDomainLongitudeValues(domainId));
    std::vector<double> vecLon {};
    for (int i = 0; i < ni; i++) {
        vecLon.push_back(-180 + (rank * ni * i) * 360 / ni_glo);
    }
    xios_handler.setDomainLongitudeValues(domainId, vecLon);
    REQUIRE(xios_handler.areDefinedDomainLongitudeValues(domainId));
    std::vector<double> vecLonOut = xios_handler.getDomainLongitudeValues(domainId);
    for (int i = 0; i < ni; i++) {
        REQUIRE(vecLonOut[i] == doctest::Approx(vecLon[i]));
    }
    // Local latitude values
    REQUIRE_FALSE(xios_handler.areDefinedDomainLatitudeValues(domainId));
    std::vector<double> vecLat {};
    for (int j = 0; j < nj; j++) {
        vecLat.push_back(-90 + j * 180 / nj_glo);
    }
    xios_handler.setDomainLatitudeValues(domainId, vecLat);
    REQUIRE(xios_handler.areDefinedDomainLatitudeValues(domainId));
    std::vector<double> vecLatOut = xios_handler.getDomainLatitudeValues(domainId);
    for (int j = 0; j < nj; j++) {
        REQUIRE(vecLatOut[j] == doctest::Approx(vecLat[j]));
    }

    // --- Tests for grid API
    std::string gridId = { "grid_2D" };
    // xios_handler.createGrid(gridId); // FIXME
    // Grid name
    REQUIRE_FALSE(xios_handler.isDefinedGridName(gridId));
    std::string gridName = { "test_grid" };
    xios_handler.setGridName(gridId, gridName);
    REQUIRE(xios_handler.isDefinedGridName(gridId));
    REQUIRE(xios_handler.getGridName(gridId) == gridName);

    // --- Tests for field API
    std::string fieldId = { "field_A" };
    // xios_handler.createField(fieldId); // FIXME
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
    // xios_handler.gridAddDomain(gridId, domainId); // FIXME and test
    // xios_handler.gridAddAxis(gridId, axisId); // FIXME and test

    // --- Tests for file API
    REQUIRE_FALSE(xios_handler.validFileId("invalid"));
    std::string fileId { "output" };
    // REQUIRE_FALSE(xios_handler.validFileId("fileId")); // TODO
    // xios_handler.createFile(fileId); // FIXME
    REQUIRE(xios_handler.validFileId(fileId));
    // File name
    REQUIRE_FALSE(xios_handler.isDefinedFileName(fileId));
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
    // xios_handler.fileAddField(fieldId, fieldId); // FIXME and test

    xios_handler.close_context_definition();

    // check the getCurrentDate method with and without ISO formatting
    xios::CDate current;
    std::string current_date = xios_handler.getCurrentDate();
    REQUIRE(current_date == "2023-03-17T17:37:00Z");
    current_date = xios_handler.getCurrentDate(false);
    REQUIRE(current_date == "2023-03-17 17:37:00");

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
        xios_handler.write(fieldId, field_A, ni, nj, axis_size);
        // verify timestep
        step = xios_handler.getCalendarStep();
        REQUIRE(step == ts);
    }

    xios_handler.context_finalize();

    // clean up
    delete[] field_A;
}
