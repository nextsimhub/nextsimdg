/*!
 * @file    XiosFile_test.cpp
 * @author  Joe Wallwork <jw2423@cam.ac.uk
 * @date    31 July 2024
 * @brief   Tests for XIOS axes
 * @details
 * This test is designed to test axis functionality of the C++ interface
 * for XIOS.
 *
 */
#include <doctest/extensions/doctest_mpi.h>
#undef INFO

#include "StructureModule/include/ParametricGrid.hpp"
#include "include/Configurator.hpp"
#include "include/Xios.hpp"

#include <iostream>

namespace Nextsim {

/*!
 * TestXiosInitialization
 *
 * This function tests the file functionality of the C++ interface for XIOS. It
 * needs to be run with 2 ranks i.e.,
 *
 * `mpirun -n 2 ./testXiosFile_MPI2`
 *
 */
MPI_TEST_CASE("TestXiosFile", 2)
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

    // Set timestep as a minimum
    xios_handler.setCalendarTimestep(Duration("P0-0T01:00:00"));

    // Axis setup
    xios_handler.createAxis("axis_A");
    xios_handler.setAxisValues("axis_A", { 0, 1 });

    // Grid setup
    xios_handler.createGrid("grid_1D");
    xios_handler.gridAddAxis("grid_1D", "axis_A");

    // Field setup
    xios_handler.createField("field_A");
    xios_handler.setFieldOperation("field_A", "instant");
    xios_handler.setFieldGridRef("field_A", "grid_1D");

    // --- Tests for file API
    const std::string fileId { "output" };
    xios_handler.createFile(fileId);
    // File name
    const std::string fileName { "diagnostic" };
    xios_handler.setFileName(fileId, fileName);
    REQUIRE(xios_handler.getFileName(fileId) == fileName);
    // File type
    const std::string fileType { "one_file" };
    xios_handler.setFileType(fileId, fileType);
    REQUIRE(xios_handler.getFileType(fileId) == fileType);
    // Output frequency
    const std::string outputFreq { "1ts" };
    xios_handler.setFileOutputFreq(fileId, outputFreq);
    REQUIRE(xios_handler.getFileOutputFreq(fileId) == outputFreq);
    // Split frequency
    const std::string splitFreq { "1ts" };
    xios_handler.setFileSplitFreq(fileId, splitFreq);
    REQUIRE(xios_handler.getFileSplitFreq(fileId) == splitFreq);
    // Add field
    xios_handler.fileAddField(fileId, "field_A");
    std::vector<std::string> fieldIds = xios_handler.fileGetFieldIds(fileId);
    REQUIRE(fieldIds.size() == 1);
    REQUIRE(fieldIds[0] == "field_A");

    xios_handler.close_context_definition();
    xios_handler.context_finalize();
}
}
