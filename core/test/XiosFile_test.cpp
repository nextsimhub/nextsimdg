/*!
 * @file    XiosFile_test.cpp
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @date    5 August 2024
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

using namespace doctest;

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
    Duration timestep("P0-0T01:30:00");
    xios_handler.setCalendarTimestep(timestep);

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
    const std::string fileId = "output";
    xios_handler.createFile(fileId);
    // File name
    const std::string fileName = "diagnostic";
    xios_handler.setFileName(fileId, fileName);
    REQUIRE(xios_handler.getFileName(fileId) == fileName);
    // File type
    const std::string fileType = "one_file";
    xios_handler.setFileType(fileId, fileType);
    REQUIRE(xios_handler.getFileType(fileId) == fileType);
    // Output frequency
    xios_handler.setFileOutputFreq(fileId, timestep);
    REQUIRE(xios_handler.getFileOutputFreq(fileId).seconds() == 1.5 * 60 * 60);
    // Split frequency
    xios_handler.setFileSplitFreq(fileId, timestep);
    REQUIRE(xios_handler.getFileSplitFreq(fileId).seconds() == 1.5 * 60 * 60);
    // File mode
    const std::string mode = "write";
    xios_handler.setFileMode(fileId, mode);
    REQUIRE(xios_handler.getFileMode(fileId) == mode);
    // File parallel access mode
    const std::string parAccess = "collective";
    xios_handler.setFileParAccess(fileId, parAccess);
    REQUIRE(xios_handler.getFileParAccess(fileId) == parAccess);
    // Add field
    xios_handler.fileAddField(fileId, "field_A");
    std::vector<std::string> fieldIds = xios_handler.fileGetFieldIds(fileId);
    REQUIRE(fieldIds.size() == 1);
    REQUIRE(fieldIds[0] == "field_A");

    // Create a new file for each time unit to check more thoroughly that XIOS interprets output
    // frequency and split frequency correctly.
    // (If we reused the same file then the XIOS interface would raise warnings.)
    xios_handler.createFile("year");
    xios_handler.setFileOutputFreq("year", Duration("P1-0T00:00:00"));
    xios_handler.setFileSplitFreq("year", Duration("P2-0T00:00:00"));
    REQUIRE(xios_handler.getFileOutputFreq("year").seconds() == 365 * 24 * 60 * 60);
    REQUIRE(xios_handler.getFileSplitFreq("year").seconds() == 2 * 365 * 24 * 60 * 60);
    xios_handler.createFile("day");
    xios_handler.setFileOutputFreq("day", Duration("P0-1T00:00:00"));
    xios_handler.setFileSplitFreq("day", Duration("P0-2T00:00:00"));
    REQUIRE(xios_handler.getFileOutputFreq("day").seconds() == 24 * 60 * 60);
    REQUIRE(xios_handler.getFileSplitFreq("day").seconds() == 2 * 24 * 60 * 60);
    xios_handler.createFile("hour");
    xios_handler.setFileOutputFreq("hour", Duration("P0-0T01:00:00"));
    xios_handler.setFileSplitFreq("hour", Duration("P0-0T02:00:00"));
    REQUIRE(xios_handler.getFileOutputFreq("hour").seconds() == 60 * 60);
    REQUIRE(xios_handler.getFileSplitFreq("hour").seconds() == 2 * 60 * 60);
    xios_handler.createFile("minute");
    xios_handler.setFileOutputFreq("minute", Duration("P0-0T00:01:00"));
    xios_handler.setFileSplitFreq("minute", Duration("P0-0T00:02:00"));
    REQUIRE(xios_handler.getFileOutputFreq("minute").seconds() == 60);
    REQUIRE(xios_handler.getFileSplitFreq("minute").seconds() == 2 * 60);
    xios_handler.createFile("second");
    xios_handler.setFileOutputFreq("second", Duration("P0-0T00:00:01"));
    xios_handler.setFileSplitFreq("second", Duration("P0-0T00:00:02"));
    REQUIRE(xios_handler.getFileOutputFreq("second").seconds() == 1);
    REQUIRE(xios_handler.getFileSplitFreq("second").seconds() == 2);

    xios_handler.close_context_definition();
    xios_handler.context_finalize();
}
}
