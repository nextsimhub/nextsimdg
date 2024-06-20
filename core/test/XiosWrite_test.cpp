/*!
 * @file    XiosWrite_test.cpp
 * @author  Joe Wallwork <jw2423@cam.ac.uk
 * @date    20 June 2024
 * @brief   Tests for XIOS write method
 * @details
 * This test is designed to test the write method of the C++ interface
 * for XIOS.
 *
 */
#include <doctest/extensions/doctest_mpi.h>

#include "include/Configurator.hpp"
#include "include/Xios.hpp"

#include <iostream>

/*!
 * TestXiosWrite
 *
 * This function tests the file writing functionality of the C++ interface for XIOS. It
 * needs to be run with 2 ranks i.e.,
 *
 * `mpirun -n 2 ./testXiosWrite_MPI2`
 *
 */
MPI_TEST_CASE("TestXiosWrite", 2)
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
    const size_t size = xios_handler.getClientMPISize();
    REQUIRE(size == 2);
    const size_t rank = xios_handler.getClientMPIRank();

    // Calendar setup
    xios_handler.setCalendarOrigin(Nextsim::TimePoint("2020-01-23T00:08:15Z"));
    xios_handler.setCalendarStart(Nextsim::TimePoint("2023-03-17T17:11:00Z"));
    xios_handler.setCalendarTimestep(Nextsim::Duration("P0-0T01:30:00"));

    // Axis setup
    xios_handler.createAxis("axis_A");
    const size_t axis_size = 30;
    xios_handler.setAxisSize("axis_A", axis_size);
    std::vector<double> axisValues(axis_size);
    for (size_t i = 0; i < axis_size; i++) {
        axisValues[i] = i;
    }
    xios_handler.setAxisValues("axis_A", axisValues);

    // Domain setup
    xios_handler.createDomain("domain_A");
    xios_handler.setDomainType("domain_A", "rectilinear");
    const size_t ni_glo = 60;
    xios_handler.setDomainGlobalLongitudeSize("domain_A", ni_glo);
    const size_t nj_glo = 20;
    xios_handler.setDomainGlobalLatitudeSize("domain_A", nj_glo);
    const size_t ni = ni_glo / size;
    xios_handler.setDomainLongitudeSize("domain_A", ni);
    const size_t nj = nj_glo;
    xios_handler.setDomainLatitudeSize("domain_A", nj);
    const size_t startLon = ni * rank;
    xios_handler.setDomainLongitudeStart("domain_A", startLon);
    const size_t startLat = 0;
    xios_handler.setDomainLatitudeStart("domain_A", startLat);
    std::vector<double> vecLon(ni);
    for (size_t i = 0; i < ni; i++) {
        vecLon[i] = -180 + (rank * ni * i) * 360 / ni_glo;
    }
    xios_handler.setDomainLongitudeValues("domain_A", vecLon);
    std::vector<double> vecLat(nj);
    for (size_t j = 0; j < nj; j++) {
        vecLat[j] = -90 + j * 180 / nj_glo;
    }
    xios_handler.setDomainLatitudeValues("domain_A", vecLat);

    // Grid setup
    xios_handler.createGrid("grid_2D");
    xios_handler.setGridName("grid_2D", "test_grid");
    xios_handler.gridAddDomain("grid_2D", "domain_A");
    xios_handler.gridAddAxis("grid_2D", "axis_A");

    // Field setup
    xios_handler.createField("field_A");
    xios_handler.setFieldName("field_A", "test_field");
    xios_handler.setFieldOperation("field_A", "instant");
    xios_handler.setFieldGridRef("field_A", "grid_2D");

    // File setup
    xios_handler.createFile("output");
    xios_handler.setFileName("output", "diagnostic");
    xios_handler.setFileType("output", "one_file");
    xios_handler.setFileOutputFreq("output", "1ts");
    xios_handler.fileAddField("output", "field_A");

    xios_handler.close_context_definition();

    // --- Tests for file API
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
        xios_handler.write("field_A", field_A, ni, nj, axis_size);
        // Verify timestep
        REQUIRE(xios_handler.getCalendarStep() == ts);
    }
    // clean up
    delete[] field_A;

    xios_handler.context_finalize();
}
