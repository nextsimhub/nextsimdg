/*!
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @brief   Tests for checking the contents of diagnostic files.
 * @details Diagnostics files should apply time-averaging.
 */
#include <doctest/extensions/doctest_mpi.h>
#undef INFO

#include "StructureModule/include/ParametricGrid.hpp"
#include "include/Finalizer.hpp"
#include "include/Model.hpp"
#include "include/ModelMPI.hpp"
#include "include/NextsimModule.hpp"
#include "include/ParaGridIO.hpp"
#include "include/Xios.hpp"
#include "include/gridNames.hpp"

#include <filesystem>

const std::string testFilesDir = TEST_FILES_DIR;
const std::string diagnosticFilename
    = testFilesDir + "/xios_test_diagnostic_20230317171100-20230317201059.nc";

static const int DGCOMP = 6;
static const int DGSTRESSCOMP = 8;
static const int CGDEGREE = 2;

namespace Nextsim {

/*!
 * TestXiosReadDiagnostic
 *
 * 1. Test the file reading functionality via `getModelState`.
 * 2. Test the file writing functionality via `writeDiagnosticTime` in the sense of checking that
 *    time-averaging is applied.
 */
MPI_TEST_CASE("TestXiosReadDiagnostic", 2)
{
    std::stringstream config;
    config << "[model]" << std::endl;
    config << "start = 2023-03-17T17:11:00Z" << std::endl;
    config << "stop = 2023-03-17T23:11:00Z" << std::endl;
    config << "time_step = P0-0T01:30:00" << std::endl;
    config << "init_file = " << diagnosticFilename << std::endl;
    config << "partition_file = xios_test_partition_metadata_2.nc" << std::endl;
    config << "[XiosInput]" << std::endl;
    config << "field_names = " << hsnowName << std::endl;
    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    // Create ModelMPI instance based off the test communicator
    auto& modelMPI = ModelMPI::getInstance(test_comm);

    // Create a Model and configure it so that time options are parsed
    // TODO: Use Model.configure for consistency with the rest of the model
    Model model;
    model.configureRestarts();
    model.configureTime();

    // Get the Xios singleton instance and check it's initialized
    // NOTE: The singleton is created when Xios::getInstance() is first called. In this test, this
    //       happens when the time sets set by ModelMetadata::setTime(). This occurs in the call to
    //       Model::configureTime() above.
    Xios& xiosHandler = Xios::getInstance();

    // Create ParametricGrid and ParaGridIO instances
    // NOTE: XIOS axes, domains, and grids are created by the ParaGridIO constructor
    Module::setImplementation<IStructure>("Nextsim::ParametricGrid");
    ParametricGrid grid;
    ParaGridIO* pio = new ParaGridIO(grid);
    grid.setIO(pio);

    xiosHandler.close_context_definition();

    // Check the input file exists
    if (!std::filesystem::exists(diagnosticFilename)) {
        throw std::runtime_error(
            "XiosReadDiagnostic_test: Input file not found. Did you run XiosWrite_test?");
    }

    // Check calendar step is zero initially
    REQUIRE(xiosHandler.getCalendarStep() == 0);

    // Deduce the local lengths of the two dimensions
    const size_t nx = ModelArray::size(ModelArray::Dimension::X);
    const size_t ny = ModelArray::size(ModelArray::Dimension::Y);

    // Read restarts from file and check they take the expected values
    for (const auto [fieldName, modelarray] : grid.getModelState(diagnosticFilename).data) {
        REQUIRE(fieldName == hsnowName);
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {
                REQUIRE(modelarray(i, j) == doctest::Approx(0.15));
            }
        }
    }

    // Remove the diagnostic files
    int rank;
    MPI_Comm_rank(test_comm, &rank);
    if (rank == 0) {
        std::filesystem::remove("xios_test_diagnostic_20230317171100-20230317201059.nc");
        std::filesystem::remove("xios_test_diagnostic_20230317201100-20230317231059.nc");
    }

    xiosHandler.context_finalize();
    Finalizer::finalize();
}
}
