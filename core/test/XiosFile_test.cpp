/*!
 * @file    XiosFile_test.cpp
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @author  Adeleke Bankole <ab3191@cam.ac.uk>
 * @date    12 May 2025
 * @brief   Tests for XIOS file
 * @details
 * This test is designed to test file functionality of the C++ interface
 * for XIOS.
 *
 */
#include <doctest/extensions/doctest_mpi.h>
#undef INFO

#include "StructureModule/include/ParametricGrid.hpp"
#include "include/Configurator.hpp"
#include "include/Finalizer.hpp"
#include "include/Xios.hpp"

using namespace doctest;

namespace Nextsim {

/*!
 * TestXiosFile
 *
 * This function tests the file functionality of the C++ interface for XIOS. It
 * needs to be run with 2 ranks i.e.,
 *
 * `mpirun -n 2 ./testXiosFile_MPI2`
 *
 */
MPI_TEST_CASE("TestXiosFile", 2)
{
    // Enable XIOS in the 'config' and provide parameters to configure it
    enableXios();
    std::stringstream config;
    config << "[model]" << std::endl;
    config << "start = 2023-03-17T17:11:00Z" << std::endl;
    config << "time_step = P0-0T01:30:00" << std::endl;
    config << "[XiosOutput]" << std::endl;
    config << "period = P0-0T03:00:00" << std::endl;
    config << "filename = output" << std::endl;
    config << "field_names = field_A" << std::endl;
    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    // Get the Xios singleton instance and check it's initialized
    Xios& xiosHandler = Xios::getInstance();
    REQUIRE(xiosHandler.isInitialized());
    const size_t size = xiosHandler.getClientMPISize();
    REQUIRE(size == 2);
    const size_t rank = xiosHandler.getClientMPIRank();

    // Create a simple axis with two points
    xiosHandler.createAxis("axis_A");
    xiosHandler.setAxisValues("axis_A", { 0.0, 1.0 });

    // Create a 1D grid comprised of the single axis
    xiosHandler.createGrid("grid_1D");
    xiosHandler.gridAddAxis("grid_1D", "axis_A");

    // Create a field on the 1D grid
    xiosHandler.createField("field_A");
    xiosHandler.setFieldOperation("field_A", "instant");
    xiosHandler.setFieldGridRef("field_A", "grid_1D");

    // --- Tests for file API
    const std::string fileId = "output";
    REQUIRE_THROWS_WITH(
        xiosHandler.getFileName("unittest_undef"), "Xios: Undefined file 'unittest_undef'");
    // File creation
    // NOTE: This is called based on the XiosInput.filename and XiosOutput.filename entries upon
    // initialization
    REQUIRE_THROWS_WITH(xiosHandler.createFile(fileId), "Xios: File 'output' already exists");
    // File name
    // NOTE: This is set based off the XiosInput.filename and XiosOutput.filename entries when a
    // file is created
    REQUIRE(xiosHandler.getFileName(fileId) == fileId);
    // File type
    REQUIRE_THROWS_WITH(xiosHandler.getFileType(fileId), "Xios: Undefined type for file 'output'");
    const std::string fileType = "one_file";
    xiosHandler.setFileType(fileId, fileType);
    REQUIRE(xiosHandler.getFileType(fileId) == fileType);
    // Output frequency
    // NOTE: This is set based off the XiosInput.period and XiosOutput.period entries when a file
    // is created
    REQUIRE(xiosHandler.getFileOutputFreq(fileId).seconds() == 3.0 * 60 * 60);
    // Split frequency
    REQUIRE_THROWS_WITH(
        xiosHandler.getFileSplitFreq(fileId), "Xios: Undefined split frequency for file 'output'");
    xiosHandler.setFileSplitFreq(fileId, xiosHandler.getCalendarTimestep());
    REQUIRE(xiosHandler.getFileSplitFreq(fileId).seconds() == 1.5 * 60 * 60);
    // File mode
    // NOTE: This is set based off the XiosInput.filename and XiosOutput.filename entries when a
    // file is created
    REQUIRE(xiosHandler.getFileMode(fileId) == "write");
    // File parallel access mode
    const std::string parAccess = "collective";
    xiosHandler.setFileParAccess(fileId, parAccess);
    REQUIRE(xiosHandler.getFileParAccess(fileId) == parAccess);
    // Check a field can be added
    xiosHandler.fileAddField(fileId, "field_A");
    std::vector<std::string> fieldIds = xiosHandler.fileGetFieldIds(fileId);
    REQUIRE(fieldIds.size() == 1);
    REQUIRE(fieldIds[0] == "field_A");

    // Create a new file for each time unit to check more thoroughly that XIOS interprets output
    // frequency and split frequency correctly.
    // (If we reused the same file then the XIOS interface would raise warnings.)
    const std::string prefix = "unittest";
    const std::string yearId = prefix + "_year";
    const std::string dayId = prefix + "_day";
    const std::string hourId = prefix + "_hour";
    const std::string minuteId = prefix + "_minute";
    const std::string secondId = prefix + "_second";
    xiosHandler.createFile(yearId);
    xiosHandler.setFileOutputFreq(yearId, Duration("P1-0T00:00:00"));
    xiosHandler.setFileSplitFreq(yearId, Duration("P2-0T00:00:00"));
    REQUIRE(xiosHandler.getFileOutputFreq(yearId).seconds() == 365 * 24 * 60 * 60);
    REQUIRE(xiosHandler.getFileSplitFreq(yearId).seconds() == 2 * 365 * 24 * 60 * 60);
    xiosHandler.createFile(dayId);
    xiosHandler.setFileOutputFreq(dayId, Duration("P0-1T00:00:00"));
    xiosHandler.setFileSplitFreq(dayId, Duration("P0-2T00:00:00"));
    REQUIRE(xiosHandler.getFileOutputFreq(dayId).seconds() == 24 * 60 * 60);
    REQUIRE(xiosHandler.getFileSplitFreq(dayId).seconds() == 2 * 24 * 60 * 60);
    xiosHandler.createFile(hourId);
    xiosHandler.setFileOutputFreq(hourId, Duration("P0-0T01:00:00"));
    xiosHandler.setFileSplitFreq(hourId, Duration("P0-0T02:00:00"));
    REQUIRE(xiosHandler.getFileOutputFreq(hourId).seconds() == 60 * 60);
    REQUIRE(xiosHandler.getFileSplitFreq(hourId).seconds() == 2 * 60 * 60);
    xiosHandler.createFile(minuteId);
    xiosHandler.setFileOutputFreq(minuteId, Duration("P0-0T00:01:00"));
    xiosHandler.setFileSplitFreq(minuteId, Duration("P0-0T00:02:00"));
    REQUIRE(xiosHandler.getFileOutputFreq(minuteId).seconds() == 60);
    REQUIRE(xiosHandler.getFileSplitFreq(minuteId).seconds() == 2 * 60);
    xiosHandler.createFile(secondId);
    xiosHandler.setFileOutputFreq(secondId, Duration("P0-0T00:00:01"));
    xiosHandler.setFileSplitFreq(secondId, Duration("P0-0T00:00:02"));
    REQUIRE(xiosHandler.getFileOutputFreq(secondId).seconds() == 1);
    REQUIRE(xiosHandler.getFileSplitFreq(secondId).seconds() == 2);

    xiosHandler.close_context_definition();
    xiosHandler.context_finalize();
    Finalizer::finalize();
}
}
