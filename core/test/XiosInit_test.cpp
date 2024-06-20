/*!
 * @file    XiosInit_test.cpp
 * @author  Joe Wallwork <jw2423@cam.ac.uk
 * @date    17 June 2024
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
    const size_t size = xios_handler.getClientMPISize();
    REQUIRE(size == 2);
    const size_t rank = xios_handler.getClientMPIRank();

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
    const std::string domainId = { "domain_A" };
    xios_handler.createDomain(domainId);
    // Domain type
    REQUIRE_FALSE(xios_handler.isDefinedDomainType(domainId));
    const std::string domainType = { "rectilinear" };
    xios_handler.setDomainType(domainId, domainType);
    REQUIRE(xios_handler.isDefinedDomainType(domainId));
    REQUIRE(xios_handler.getDomainType(domainId) == domainType);
    // Global longitude size
    REQUIRE_FALSE(xios_handler.isDefinedDomainGlobalLongitudeSize(domainId));
    const size_t ni_glo = 60;
    xios_handler.setDomainGlobalLongitudeSize(domainId, ni_glo);
    REQUIRE(xios_handler.isDefinedDomainGlobalLongitudeSize(domainId));
    REQUIRE(xios_handler.getDomainGlobalLongitudeSize(domainId) == ni_glo);
    // Global latitude size
    REQUIRE_FALSE(xios_handler.isDefinedDomainGlobalLatitudeSize(domainId));
    const size_t nj_glo = 20;
    xios_handler.setDomainGlobalLatitudeSize(domainId, nj_glo);
    REQUIRE(xios_handler.isDefinedDomainGlobalLatitudeSize(domainId));
    REQUIRE(xios_handler.getDomainGlobalLatitudeSize(domainId) == nj_glo);
    // Local longitude size
    REQUIRE_FALSE(xios_handler.isDefinedDomainLongitudeSize(domainId));
    const size_t ni = ni_glo / size;
    xios_handler.setDomainLongitudeSize(domainId, ni);
    REQUIRE_FALSE(xios_handler.isDefinedDomainLatitudeSize(domainId));
    REQUIRE(xios_handler.getDomainLongitudeSize(domainId) == ni);
    // Local latitude size
    REQUIRE_FALSE(xios_handler.isDefinedDomainLatitudeSize(domainId));
    const size_t nj = nj_glo;
    xios_handler.setDomainLatitudeSize(domainId, nj);
    REQUIRE(xios_handler.isDefinedDomainLatitudeSize(domainId));
    REQUIRE(xios_handler.getDomainLatitudeSize(domainId) == nj);
    // Local longitude start
    REQUIRE_FALSE(xios_handler.isDefinedDomainLongitudeStart(domainId));
    const size_t startLon = ni * rank;
    xios_handler.setDomainLongitudeStart(domainId, startLon);
    REQUIRE(xios_handler.isDefinedDomainLongitudeStart(domainId));
    REQUIRE(xios_handler.getDomainLongitudeStart(domainId) == startLon);
    // Local latitude start
    REQUIRE_FALSE(xios_handler.isDefinedDomainLatitudeStart(domainId));
    const size_t startLat = 0;
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
    const std::string gridId = { "grid_2D" };
    xios_handler.createGrid(gridId);
    // Grid name
    REQUIRE_FALSE(xios_handler.isDefinedGridName(gridId));
    const std::string gridName = { "test_grid" };
    xios_handler.setGridName(gridId, gridName);
    REQUIRE(xios_handler.isDefinedGridName(gridId));
    REQUIRE(xios_handler.getGridName(gridId) == gridName);

    // --- Tests for field API
    const std::string fieldId = { "field_A" };
    xios_handler.createField(fieldId);
    // Field name
    const std::string fieldName = { "test_field" };
    REQUIRE_FALSE(xios_handler.isDefinedFieldName(fieldId));
    xios_handler.setFieldName(fieldId, fieldName);
    REQUIRE(xios_handler.getFieldName(fieldId) == fieldName);
    REQUIRE(xios_handler.isDefinedFieldName(fieldId));
    // Operation
    const std::string operation = { "instant" };
    REQUIRE_FALSE(xios_handler.isDefinedFieldOperation(fieldId));
    xios_handler.setFieldOperation(fieldId, operation);
    REQUIRE(xios_handler.isDefinedFieldOperation(fieldId));
    REQUIRE(xios_handler.getFieldOperation(fieldId) == operation);
    // Grid reference
    const std::string gridRef = { "grid_2D" };
    REQUIRE_FALSE(xios_handler.isDefinedFieldGridRef(fieldId));
    xios_handler.setFieldGridRef(fieldId, gridRef);
    REQUIRE(xios_handler.isDefinedFieldGridRef(fieldId));
    REQUIRE(xios_handler.getFieldGridRef(fieldId) == gridRef);
    xios_handler.gridAddDomain(gridId, domainId);
    xios_handler.gridAddAxis(gridId, axisId);

    // --- Tests for file API
    const std::string fileId { "output" };
    REQUIRE_FALSE(xios_handler.validFileId(fileId));
    xios_handler.createFile(fileId);
    REQUIRE(xios_handler.validFileId(fileId));
    // File name
    const std::string fileName { "diagnostic" };
    xios_handler.setFileName(fileId, fileName);
    REQUIRE(xios_handler.isDefinedFileName(fileId));
    REQUIRE(xios_handler.getFileName(fileId) == fileName);
    // File type
    const std::string fileType { "one_file" };
    REQUIRE_FALSE(xios_handler.isDefinedFileType(fileId));
    xios_handler.setFileType(fileId, fileType);
    REQUIRE(xios_handler.isDefinedFileType(fileId));
    REQUIRE(xios_handler.getFileType(fileId) == fileType);
    // Output frequency
    REQUIRE_FALSE(xios_handler.isDefinedFileOutputFreq(fileId));
    const std::string freq { "1ts" };
    xios_handler.setFileOutputFreq(fileId, freq);
    REQUIRE(xios_handler.isDefinedFileOutputFreq(fileId));
    REQUIRE(xios_handler.getFileOutputFreq(fileId) == freq);
    xios_handler.fileAddField(fileId, fieldId);

    xios_handler.close_context_definition();

    // create some fake data to test writing methods
    double* field_A = new double[ni * nj * axis_size];
    for (size_t idx = 0; idx < ni * nj * axis_size; idx++) {
        field_A[idx] = 1.0 * idx;
    }

    // Verify calendar step is starting from zero
    REQUIRE(xios_handler.getCalendarStep() == 0);

    // simulate 4 iterations (timesteps)
    for (int ts = 1; ts <= 4; ts++) {
        // update the current timestep
        xios_handler.updateCalendar(ts);
        // send data to XIOS to be written to disk
        xios_handler.write(fieldId, field_A, ni, nj, axis_size);
        // Verify timestep
        REQUIRE(xios_handler.getCalendarStep() == ts);
    }

    xios_handler.context_finalize();

    // clean up
    delete[] field_A;
}
