/*!
 * @file    XiosFile_test.cpp
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @date    02 Jan 2025
 * @brief   Tests for XIOS axes
 * @details
 * This test is designed to test axis functionality of the C++ interface
 * for XIOS.
 *
 */
#include <doctest/extensions/doctest_mpi.h>
#undef INFO

#include "StructureModule/include/ParametricGrid.hpp"
#include "include/Xios.hpp"

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
    enableXios();

    // Get the Xios singleton instance and check it's initialized
    Xios* xiosHandler = Xios::getInstance("P0-0T01:30:00");
    REQUIRE(xiosHandler->isInitialized());
    const size_t size = xiosHandler->getClientMPISize();
    REQUIRE(size == 2);
    const size_t rank = xiosHandler->getClientMPIRank();

    // Create a simple axis with two points
    xiosHandler->createAxis("axis_A");
    xiosHandler->setAxisValues("axis_A", { 0.0, 1.0 });

    // Create a 1D grid comprised of the single axis
    xiosHandler->createGrid("grid_1D");
    xiosHandler->gridAddAxis("grid_1D", "axis_A");

    // Create a field on the 1D grid
    xiosHandler->createField("field_A");
    xiosHandler->setFieldOperation("field_A", "instant");
    xiosHandler->setFieldGridRef("field_A", "grid_1D");

    // --- Tests for file API
    const std::string fileId = "output";
    REQUIRE_THROWS_WITH(xiosHandler->getFileName(fileId), "Xios: Undefined file 'output'");
    xiosHandler->createFile(fileId);
    REQUIRE_THROWS_WITH(xiosHandler->createFile(fileId), "Xios: File 'output' already exists");
    // File name
    REQUIRE_THROWS_WITH(xiosHandler->getFileName(fileId), "Xios: Undefined name for file 'output'");
    const std::string fileName = "diagnostic";
    xiosHandler->setFileName(fileId, fileName);
    REQUIRE(xiosHandler->getFileName(fileId) == fileName);
    // File type
    REQUIRE_THROWS_WITH(xiosHandler->getFileType(fileId), "Xios: Undefined type for file 'output'");
    const std::string fileType = "one_file";
    xiosHandler->setFileType(fileId, fileType);
    REQUIRE(xiosHandler->getFileType(fileId) == fileType);
    // Output frequency
    REQUIRE_THROWS_WITH(xiosHandler->getFileOutputFreq(fileId),
        "Xios: Undefined output frequency for file 'output'");
    Duration timestep = xiosHandler->getCalendarTimestep();
    xiosHandler->setFileOutputFreq(fileId, timestep);
    REQUIRE(xiosHandler->getFileOutputFreq(fileId).seconds() == 1.5 * 60 * 60);
    // Split frequency
    REQUIRE_THROWS_WITH(
        xiosHandler->getFileSplitFreq(fileId), "Xios: Undefined split frequency for file 'output'");
    xiosHandler->setFileSplitFreq(fileId, timestep);
    REQUIRE(xiosHandler->getFileSplitFreq(fileId).seconds() == 1.5 * 60 * 60);
    // File mode
    const std::string mode = "write";
    xiosHandler->setFileMode(fileId, mode);
    REQUIRE(xiosHandler->getFileMode(fileId) == mode);
    // File parallel access mode
    const std::string parAccess = "collective";
    xiosHandler->setFileParAccess(fileId, parAccess);
    REQUIRE(xiosHandler->getFileParAccess(fileId) == parAccess);
    // Add field
    xiosHandler->fileAddField(fileId, "field_A");
    std::vector<std::string> fieldIds = xiosHandler->fileGetFieldIds(fileId);
    REQUIRE(fieldIds.size() == 1);
    REQUIRE(fieldIds[0] == "field_A");

    // Create a new file for each time unit to check more thoroughly that XIOS interprets output
    // frequency and split frequency correctly.
    // (If we reused the same file then the XIOS interface would raise warnings.)
    xiosHandler->createFile("year");
    xiosHandler->setFileOutputFreq("year", Duration("P1-0T00:00:00"));
    xiosHandler->setFileSplitFreq("year", Duration("P2-0T00:00:00"));
    REQUIRE(xiosHandler->getFileOutputFreq("year").seconds() == 365 * 24 * 60 * 60);
    REQUIRE(xiosHandler->getFileSplitFreq("year").seconds() == 2 * 365 * 24 * 60 * 60);
    xiosHandler->createFile("day");
    xiosHandler->setFileOutputFreq("day", Duration("P0-1T00:00:00"));
    xiosHandler->setFileSplitFreq("day", Duration("P0-2T00:00:00"));
    REQUIRE(xiosHandler->getFileOutputFreq("day").seconds() == 24 * 60 * 60);
    REQUIRE(xiosHandler->getFileSplitFreq("day").seconds() == 2 * 24 * 60 * 60);
    xiosHandler->createFile("hour");
    xiosHandler->setFileOutputFreq("hour", Duration("P0-0T01:00:00"));
    xiosHandler->setFileSplitFreq("hour", Duration("P0-0T02:00:00"));
    REQUIRE(xiosHandler->getFileOutputFreq("hour").seconds() == 60 * 60);
    REQUIRE(xiosHandler->getFileSplitFreq("hour").seconds() == 2 * 60 * 60);
    xiosHandler->createFile("minute");
    xiosHandler->setFileOutputFreq("minute", Duration("P0-0T00:01:00"));
    xiosHandler->setFileSplitFreq("minute", Duration("P0-0T00:02:00"));
    REQUIRE(xiosHandler->getFileOutputFreq("minute").seconds() == 60);
    REQUIRE(xiosHandler->getFileSplitFreq("minute").seconds() == 2 * 60);
    xiosHandler->createFile("second");
    xiosHandler->setFileOutputFreq("second", Duration("P0-0T00:00:01"));
    xiosHandler->setFileSplitFreq("second", Duration("P0-0T00:00:02"));
    REQUIRE(xiosHandler->getFileOutputFreq("second").seconds() == 1);
    REQUIRE(xiosHandler->getFileSplitFreq("second").seconds() == 2);

    xiosHandler->close_context_definition();
    xiosHandler->context_finalize();
}
}
