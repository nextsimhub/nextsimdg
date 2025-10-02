/*!
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @author  Adeleke Bankole <ab3191@cam.ac.uk>
 * @brief   Tests for XIOS file
 * @details
 * This test is designed to test file functionality of the C++ interface
 * for XIOS.
 */
#include <doctest/extensions/doctest_mpi.h>
#undef INFO

#include "StructureModule/include/ParametricGrid.hpp"
#include "include/Configurator.hpp"
#include "include/Finalizer.hpp"
#include "include/ModelMPI.hpp"
#include "include/ModelMetadata.hpp"
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
    config << "[XiosInput]" << std::endl;
    config << "period = P0-0T03:00:00" << std::endl;
    config << "filename = xios_test_input.nc" << std::endl;
    config << "field_names = field_1D" << std::endl;
    config << "[XiosOutput]" << std::endl;
    config << "period = P0-0T03:00:00" << std::endl;
    config << "filename = xios_test_output.nc" << std::endl;
    config << "field_names = field_2D" << std::endl;
    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    // Create ModelMetadata instance based off a partition metadata file
    auto& modelMPI = ModelMPI::getInstance(test_comm);
    auto& metadata = ModelMetadata::getInstance("xios_test_partition_metadata_2.nc");

    // Get the Xios singleton instance and check it's initialized
    Xios& xiosHandler = Xios::getInstance();
    xiosHandler.affixModelMetadata();
    REQUIRE(xiosHandler.isInitialized());
    REQUIRE(xiosHandler.getClientMPISize() == 2);

    // Create a vertical axis, too
    xiosHandler.createAxis("z_axis");
    xiosHandler.setAxisValues("z_axis", { 0.0, 1.0 });

    // Create a 1D grid
    // NOTE: The 2D grid is created automatically along with the xy_domain
    xiosHandler.createGrid("grid_1D");
    xiosHandler.gridAddAxis("grid_1D", "z_axis");

    // Associate fields with grids
    // NOTE: fields are automatically created along with files
    xiosHandler.setFieldOperation("field_1D", "instant");
    xiosHandler.setFieldGridRef("field_1D", "grid_1D");
    xiosHandler.setFieldOperation("field_2D", "instant");
    xiosHandler.setFieldGridRef("field_2D", "grid_2D");

    // --- Tests for file API
    const std::string inFileId = "xios_test_input";
    const std::string outFileId = "xios_test_output";
    REQUIRE_THROWS_WITH(
        xiosHandler.getFileType("unittest_undef"), "Xios: Undefined file 'unittest_undef'");
    // File creation
    // NOTE: This is called based on the XiosInput.filename and XiosOutput.filename entries upon
    // initialization
    REQUIRE_THROWS_WITH(
        xiosHandler.createFile(outFileId), "Xios: File 'xios_test_output' already exists");
    // File type
    // NOTE: This is to "one_file" when createFile is called
    REQUIRE(xiosHandler.getFileType(outFileId) == "one_file");
    // Output frequency
    // NOTE: This is set based off the XiosInput.period and XiosOutput.period entries when a file
    // is created
    REQUIRE(xiosHandler.getFileOutputFreq(outFileId).seconds() == 3.0 * 60 * 60);
    // Split frequency
    REQUIRE_THROWS_WITH(xiosHandler.getFileSplitFreq(outFileId),
        "Xios: Undefined split frequency for file 'xios_test_output'");
    xiosHandler.setFileSplitFreq(outFileId, xiosHandler.getCalendarTimestep());
    REQUIRE(xiosHandler.getFileSplitFreq(outFileId).seconds() == 1.5 * 60 * 60);
    // File mode
    // NOTE: setFileMode is set based off the XiosInput.filename and XiosOutput.filename entries
    // when a file is created
    REQUIRE(xiosHandler.getFileMode(inFileId) == "read");
    REQUIRE(xiosHandler.getFileMode(outFileId) == "write");
    // File parallel access mode
    // NOTE: setFileParAccess is is to "collective" when a file is created for reading
    REQUIRE(xiosHandler.getFileParAccess(inFileId) == "collective");
    // File add field
    // NOTE: fileAddField is triggered by a call to createFile, which parses the config to create
    // all the corresponding fields and then associated them with the file
    std::vector<std::string> field_1DIds = xiosHandler.fileGetFieldIds(inFileId);
    REQUIRE(field_1DIds.size() == 1);
    REQUIRE(field_1DIds[0] == "field_1D");
    std::vector<std::string> field_2DIds = xiosHandler.fileGetFieldIds(outFileId);
    REQUIRE(field_2DIds.size() == 1);
    REQUIRE(field_2DIds[0] == "field_2D");

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
