/*!
 * @file    XiosWrite_test.cpp
 * @author  Joe Wallwork <jw2423@cam.ac.uk
 * @date    27 June 2024
 * @brief   Tests for XIOS write method
 * @details
 * This test is designed to test the write method of the C++ interface
 * for XIOS.
 *
 */
#include <doctest/extensions/doctest_mpi.h>
#undef INFO

#include "include/Configurator.hpp"
// #include "include/ConfiguredModule.hpp"
#include "include/ParametricGrid.hpp"
#include "include/StructuredModule.hpp"
#include "include/Xios.hpp"

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
    xios_handler.createFile("output");
    xios_handler.setFileType("output", "one_file");
    xios_handler.setFileOutputFreq("output", "1ts");
    xios_handler.fileAddField("output", "field_2D");
    xios_handler.fileAddField("output", "field_3D");
    xios_handler.fileAddField("output", "field_4D");

    xios_handler.close_context_definition();

    // --- Tests for writing to file
    Module::setImplementation<IStructure>("Nextsim::ParametricGrid");
    ParametricGrid grid;
    ModelArray::setDimension(ModelArray::Dimension::X, ni);
    ModelArray::setDimension(ModelArray::Dimension::Y, nj);
    ModelArray::setDimension(ModelArray::Dimension::Z, axis_size);
    // Create some fake data to test writing methods
    HField field_A(ModelArray::Type::Z);
    field_A.resize();
    for (size_t k = 0; k < axis_size; ++k) {
        for (size_t j = 0; j < nj; ++j) {
            for (size_t i = 0; i < ni; ++i) {
                field_A(i, j, k) = 1.0 * (i + ni * (j + nj * k));
            }
        }
    }
    // Verify calendar step is starting from zero
    int step = xios_handler.getCalendarStep();
    REQUIRE(step == 0);
    // Simulate 4 iterations (timesteps)
    for (int ts = 1; ts <= 4; ts++) {
        // Update the current timestep
        xios_handler.updateCalendar(ts);
        // Send data to XIOS to be written to disk
        xios_handler.write(fieldId, field_A);
        // Verify timestep
        step = xios_handler.getCalendarStep();
        REQUIRE(step == ts);
    }

    xios_handler.context_finalize();
}
}
