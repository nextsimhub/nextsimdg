/*!
 * @file    XiosWrite_test.cpp
 * @author  Joe Wallwork <jw2423@cam.ac.uk
 * @date    31 July 2024
 * @brief   Tests for XIOS write method
 * @details
 * This test is designed to test the write method of the C++ interface
 * for XIOS.
 *
 */
#include <doctest/extensions/doctest_mpi.h>
#undef INFO

#include "StructureModule/include/ParametricGrid.hpp"
#include "include/Configurator.hpp"
#include "include/Module.hpp"
#include "include/Xios.hpp"

#include <filesystem>
#include <iostream>

namespace Nextsim {

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
    xios_handler.createFile("xios_test_output");
    xios_handler.setFileType("xios_test_output", "one_file");
    xios_handler.setFileOutputFreq("xios_test_output", "1ts");
    xios_handler.setFileSplitFreq("xios_test_output", "2ts");
    xios_handler.fileAddField("xios_test_output", "field_2D");
    xios_handler.fileAddField("xios_test_output", "field_3D");
    xios_handler.fileAddField("xios_test_output", "field_4D");

    xios_handler.close_context_definition();

    // --- Tests for writing to file
    Module::setImplementation<IStructure>("Nextsim::ParametricGrid");
    ModelArray::setDimension(ModelArray::Dimension::X, n1);
    ModelArray::setDimension(ModelArray::Dimension::Y, n2);
    ModelArray::setDimension(ModelArray::Dimension::Z, n3);
    // Create some fake data to test writing methods
    HField field_2D(ModelArray::Type::H);
    field_2D.resize();
    for (size_t j = 0; j < n2; ++j) {
        for (size_t i = 0; i < n1; ++i) {
            field_2D(i, j) = 1.0 * (i + n1 * j);
        }
    }
    HField field_3D(ModelArray::Type::Z);
    field_3D.resize();
    for (size_t k = 0; k < n3; ++k) {
        for (size_t j = 0; j < n2; ++j) {
            for (size_t i = 0; i < n1; ++i) {
                field_3D(i, j, k) = 1.0 * (i + n1 * (j + n2 * k));
            }
        }
    }
    // TODO: field_4D?
    // Verify calendar step is starting from zero
    REQUIRE(xios_handler.getCalendarStep() == 0);
    // Check a file with the expected name doesn't exist yet
    REQUIRE_FALSE(std::filesystem::exists("xios_test_output*.nc"));
    // Simulate 4 iterations (timesteps)
    for (int ts = 1; ts <= 4; ts++) {
        // Update the current timestep
        xios_handler.updateCalendar(ts);
        // Send data to XIOS to be written to disk
        xios_handler.write("field_2D", field_2D);
        xios_handler.write("field_3D", field_3D);
        // TODO: field_4D?
        // Verify timestep
        REQUIRE(xios_handler.getCalendarStep() == ts);
    }
    // Check the files have indeed been created then remove it
    // FIXME: These aren't the datetimes that come out - seems to be using 54hr timestep
    REQUIRE(std::filesystem::exists("xios_test_output_20230317171100-20230317201059.nc"));
    REQUIRE(std::filesystem::exists("xios_test_output_20230317201100-20230317231059.nc"));
    if (rank == 0) {
        std::filesystem::remove("xios_test_output_20230317171100-20230317201059.nc");
        std::filesystem::remove("xios_test_output_20230317201100-20230317231059.nc");
    }
    xios_handler.context_finalize();
}
}
