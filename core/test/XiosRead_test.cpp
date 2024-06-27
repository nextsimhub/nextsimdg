/*!
 * @file    XiosRead_test.cpp
 * @author  Joe Wallwork <jw2423@cam.ac.uk
 * @date    27 June 2024
 * @brief   Tests for XIOS read method
 * @details
 * This test is designed to test the read method of the C++ interface
 * for XIOS.
 *
 */
#include <doctest/extensions/doctest_mpi.h>
#undef INFO

#include "include/Configurator.hpp"
#include "include/Xios.hpp"

#include <filesystem>
#include <iostream>

namespace Nextsim {

/*!
 * TestXiosRead
 *
 * This function tests the file reading functionality of the C++ interface for XIOS. It
 * needs to be run with 2 ranks i.e.,
 *
 * `mpirun -n 2 ./testXiosRead_MPI2`
 *
 */
MPI_TEST_CASE("TestXiosRead", 2)
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
    const size_t size = xios_handler.getClientMPISize();
    REQUIRE(size == 2);
    const size_t rank = xios_handler.getClientMPIRank();

    // Calendar setup
    xios_handler.setCalendarOrigin(TimePoint("2020-01-23T00:08:15Z"));
    xios_handler.setCalendarStart(TimePoint("2023-03-17T17:11:00Z"));
    xios_handler.setCalendarTimestep(Duration("P0-0T01:30:00"));

    // Axis setup
    const int n1 = 2;
    const int n2 = 3;
    const int n3 = 4;
    const int n4 = 5;
    xios_handler.createAxis("axis_A");
    xios_handler.setAxisSize("axis_A", n1);
    xios_handler.setAxisValues("axis_A", { 0, 1 });
    xios_handler.createAxis("axis_B");
    xios_handler.setAxisSize("axis_B", n2);
    xios_handler.setAxisValues("axis_B", { 0, 1, 2 });
    xios_handler.createAxis("axis_C");
    xios_handler.setAxisSize("axis_C", n3);
    xios_handler.setAxisValues("axis_C", { 0, 1, 2, 3 });
    xios_handler.createAxis("axis_D");
    xios_handler.setAxisSize("axis_D", n4);
    xios_handler.setAxisValues("axis_D", { 0, 1, 2, 3, 4 });

    // Grid setup
    xios_handler.createGrid("grid_2D");
    xios_handler.gridAddAxis("grid_2D", "axis_A");
    xios_handler.gridAddAxis("grid_2D", "axis_B");
    xios_handler.createGrid("grid_3D");
    xios_handler.gridAddAxis("grid_3D", "axis_A");
    xios_handler.gridAddAxis("grid_3D", "axis_B");
    xios_handler.gridAddAxis("grid_3D", "axis_C");
    xios_handler.createGrid("grid_4D");
    xios_handler.gridAddAxis("grid_4D", "axis_A");
    xios_handler.gridAddAxis("grid_4D", "axis_B");
    xios_handler.gridAddAxis("grid_4D", "axis_C");
    xios_handler.gridAddAxis("grid_4D", "axis_D");

    // Field setup
    xios_handler.createField("field_2D");
    xios_handler.setFieldOperation("field_2D", "instant");
    xios_handler.setFieldGridRef("field_2D", "grid_2D");
    xios_handler.createField("field_3D");
    xios_handler.setFieldOperation("field_3D", "instant");
    xios_handler.setFieldGridRef("field_3D", "grid_3D");
    xios_handler.createField("field_4D");
    xios_handler.setFieldOperation("field_4D", "instant");
    xios_handler.setFieldGridRef("field_4D", "grid_4D");

    // File setup
    xios_handler.createFile("xios_test_input");
    xios_handler.setFileType("xios_test_input", "one_file");
    xios_handler.setFileOutputFreq("xios_test_input", "1ts");
    xios_handler.setFileMode("xios_test_input", "read");
    xios_handler.setFileParAccess("xios_test_input", "collective");
    xios_handler.fileAddField("xios_test_input", "field_2D");
    xios_handler.fileAddField("xios_test_input", "field_3D");
    xios_handler.fileAddField("xios_test_input", "field_4D");

    xios_handler.close_context_definition();

    // --- Tests for file API
    // create arrays for data to be read into
    double* field_2D = new double[n1 * n2];
    double* field_3D = new double[n1 * n2 * n3];
    double* field_4D = new double[n1 * n2 * n3 * n4];
    // Verify calendar step is starting from zero
    REQUIRE(xios_handler.getCalendarStep() == 0);
    // Check the input file exists
    REQUIRE(std::filesystem::exists("xios_test_input.nc"));
    // Simulate 4 iterations (timesteps)
    for (int ts = 1; ts <= 4; ts++) {
        // update the current timestep
        xios_handler.updateCalendar(ts);
        // receive data from XIOS that is read from disk
        xios_handler.read("field_2D", field_2D, n1, n2);
        xios_handler.read("field_3D", field_3D, n1, n2, n3);
        xios_handler.read("field_4D", field_4D, n1, n2, n3, n4);
        // Verify timestep
        REQUIRE(xios_handler.getCalendarStep() == ts);
    }
    // Verify fields have been read in correctly
    for (size_t idx = 0; idx < n1 * n2; idx++) {
        REQUIRE(field_2D[idx] == doctest::Approx(idx));
    }
    for (size_t idx = 0; idx < n1 * n2 * n3; idx++) {
        REQUIRE(field_3D[idx] == doctest::Approx(idx));
    }
    for (size_t idx = 0; idx < n1 * n2 * n3 * n4; idx++) {
        REQUIRE(field_4D[idx] == doctest::Approx(idx));
    }
    // Clean up
    delete[] field_2D;
    delete[] field_3D;
    delete[] field_4D;

    xios_handler.context_finalize();
}
}
