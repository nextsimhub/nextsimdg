/*!
 * @file    XiosRead_test.cpp
 * @author  Joe Wallwork <jw2423@cam.ac.uk
 * @date    21 August 2024
 * @brief   Tests for XIOS read method
 * @details
 * This test is designed to test the read method of the C++ interface
 * for XIOS.
 *
 */
#include <doctest/extensions/doctest_mpi.h>
#undef INFO

#include "StructureModule/include/ParametricGrid.hpp"
#include "include/Configurator.hpp"
#include "include/NextsimModule.hpp"
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
    Duration timestep("P0-0T01:30:00");
    xios_handler.setCalendarTimestep(timestep);

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
    xios_handler.setFieldReadAccess("field_2D", true);
    xios_handler.setFieldFreqOffset("field_2D", timestep);
    xios_handler.createField("field_3D");
    xios_handler.setFieldOperation("field_3D", "instant");
    xios_handler.setFieldGridRef("field_3D", "grid_3D");
    xios_handler.setFieldReadAccess("field_3D", true);
    xios_handler.setFieldFreqOffset("field_3D", timestep);
    xios_handler.createField("field_4D");
    xios_handler.setFieldOperation("field_4D", "instant");
    xios_handler.setFieldGridRef("field_4D", "grid_4D");
    xios_handler.setFieldReadAccess("field_4D", true);
    xios_handler.setFieldFreqOffset("field_4D", timestep);

    // File setup
    xios_handler.createFile("xios_test_input");
    xios_handler.setFileType("xios_test_input", "one_file");
    xios_handler.setFileOutputFreq("xios_test_input", timestep);
    xios_handler.setFileMode("xios_test_input", "read");
    xios_handler.setFileParAccess("xios_test_input", "collective");
    xios_handler.fileAddField("xios_test_input", "field_2D");
    xios_handler.fileAddField("xios_test_input", "field_3D");
    xios_handler.fileAddField("xios_test_input", "field_4D");

    xios_handler.close_context_definition();

    // --- Tests for reading to file
    Module::setImplementation<IStructure>("Nextsim::ParametricGrid");
    ModelArray::setDimension(ModelArray::Dimension::X, n1, n1, 0);
    ModelArray::setDimension(ModelArray::Dimension::Y, n2, n2, 0);
    ModelArray::setDimension(ModelArray::Dimension::Z, n3, n3, 0);
    // Create some fake data to test writing methods
    HField field_2D(ModelArray::Type::H);
    field_2D.resize();
    HField field_3D(ModelArray::Type::Z);
    field_3D.resize();
    // TODO: Implement 4D case
    // Verify calendar step is starting from zero
    REQUIRE(xios_handler.getCalendarStep() == 0);
    // Check the input file exists
    REQUIRE(std::filesystem::exists("xios_test_input.nc"));
    // Simulate 4 iterations (timesteps)
    for (int ts = 1; ts <= 4; ts++) {
        // Update the current timestep
        xios_handler.updateCalendar(ts);
        // Receive data from XIOS that is read from disk
        xios_handler.read("field_2D", field_2D);
        xios_handler.read("field_3D", field_3D);
        // TODO: Implement 4D case
        // Verify timestep
        REQUIRE(xios_handler.getCalendarStep() == ts);
    }
    // Verify fields have been read in correctly
    for (size_t j = 0; j < n2; ++j) {
        for (size_t i = 0; i < n1; ++i) {
            REQUIRE(field_2D(i, j) == doctest::Approx(1.0 * (i + n1 * j)));
        }
    }
    for (size_t k = 0; k < n3; ++k) {
        for (size_t j = 0; j < n2; ++j) {
            for (size_t i = 0; i < n1; ++i) {
                REQUIRE(field_3D(i, j, k) == doctest::Approx(1.0 * (i + n1 * (j + n2 * k))));
            }
        }
    }

    xios_handler.context_finalize();
}
}
