/*!
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @brief   Tests for checking the contents of diagnostic files.
 * @details Diagnostics files should apply time-averaging.
 */
#include <doctest/extensions/doctest_mpi.h>
#undef INFO

#include "include/Finalizer.hpp"
#include "include/Model.hpp"
#include "include/ModelMPI.hpp"
#include "include/NextsimModule.hpp"
#include "include/StructureFactory.hpp"
#include "include/Xios.hpp"
#include "include/gridNames.hpp"

#include <filesystem>

const std::string testFilesDir = TEST_FILES_DIR;
const std::string diagnosticFilename
    = testFilesDir + "/diagnostic_2023-03-17T00:00:00Z-2023-03-17T02:59:59Z.nc";

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
MPI_TEST_CASE("TestXiosReadDiagnostic", 3)
{
    std::stringstream config;
    config << "[model]" << std::endl;
    config << "start = 2023-03-17T00:00:00Z" << std::endl;
    config << "stop = 2023-03-17T00:00:00Z" << std::endl;
    config << "time_step = P0-0T01:30:00" << std::endl;
    config << "init_file = " << diagnosticFilename << std::endl;
    config << "restart_period = P0-0T03:00:00" << std::endl;
    config << "partition_file = xios_test_partition_metadata_3.nc" << std::endl;
    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    // Check the input file exists
    if (!std::filesystem::exists(diagnosticFilename)) {
        throw std::runtime_error(
            "XiosReadDiagnostic_test: Input file not found. Did you run XiosWriteDiagnostic_test?");
    }

    // Create ModelMPI instance based off the test communicator
    auto& modelMPI = ModelMPI::getInstance(test_comm);

    // Create a Model and configure it so that time options are parsed
    Model model;
    model.configureRestarts();
    model.configureTime();

    // Get the Xios singleton instance and check calendar step is zero initially
    // NOTE: The singleton is created when Xios::getInstance() is first called. In this test, this
    //       happens when the time sets set by ModelMetadata::setTime(). This occurs in the call to
    //       Model::configureTime() above.
    Xios& xiosHandler = Xios::getInstance();
    REQUIRE(xiosHandler.getCalendarStep() == 0);

    // Read restarts from file and check they take the expected values
    // NOTE: The ParametricGrid is created and the XIOS context definition is closed in the call to
    //       StructureFactory::stateFromFile()
    ModelState modelstate = StructureFactory::stateFromFile(diagnosticFilename);
    REQUIRE(modelstate.data.count(sssName) > 0);
    const size_t nx = ModelArray::size(ModelArray::Dimension::X);
    const size_t ny = ModelArray::size(ModelArray::Dimension::Y);
    REQUIRE(modelstate.data.count(sssName) > 0);
    REQUIRE(modelstate.data.count(latitudeName) > 0);
    REQUIRE(modelstate.data.count(longitudeName) > 0);
    for (const auto [fieldName, modelarray] : modelstate.data) {
        if (fieldName == sssName) {
            for (size_t j = 0; j < ny; ++j) {
                for (size_t i = 0; i < nx; ++i) {
                    REQUIRE(modelarray(i, j) == doctest::Approx(0.15));
                }
            }
        }
    }

    // Remove the diagnostic files
    int rank;
    MPI_Comm_rank(test_comm, &rank);
    if (rank == 0) {
        std::filesystem::remove("diagnostic_2023-03-17T00:00:00Z-2023-03-17T02:59:59Z.nc");
        std::filesystem::remove("diagnostic_2023-03-17T03:00:00Z-2023-03-17T05:59:59Z.nc");
    }

    xiosHandler.context_finalize();
    Finalizer::finalize();
}
}
