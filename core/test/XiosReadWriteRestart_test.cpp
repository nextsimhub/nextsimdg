/*!
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @brief   Tests for XIOS functionality for reading and writing restart files.
 * @details The functionality of both reading and writing restarts via XIOS is tested.
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
const std::string inputFilename = testFilesDir + "/xios_test_input.nc";
const std::string outputFilename = testFilesDir + "/readwrite%Y-%m-%dT%H:%M:%SZ.nc";

static const int DGCOMP = 6;
static const int DGSTRESSCOMP = 8;
static const int CGDEGREE = 2;

namespace Nextsim {

/*!
 * TestXiosReadWriteRestart
 *
 * 1. Test reading of restarts via `getModelState`.
 * 2. Test writing of restarts via `dumpModelState`.
 */
MPI_TEST_CASE("TestXiosReadWriteRestart", 2)
{
    std::stringstream config;
    config << "[model]" << std::endl;
    config << "start = 2023-03-17T17:11:00Z" << std::endl;
    config << "stop = 2023-03-17T23:11:00Z" << std::endl;
    config << "time_step = P0-0T01:30:00" << std::endl;
    config << "init_file = " << inputFilename << std::endl;
    config << "restart_file = " << outputFilename << std::endl;
    config << "partition_file = xios_test_partition_metadata_2.nc" << std::endl;
    config << "restart_period = P0-0T03:00:00" << std::endl;
    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    // Create ModelMPI instance based off the test communicator
    auto& modelMPI = ModelMPI::getInstance(test_comm);

    // Create a Model and configure it so that time options are parsed
    Model model;
    model.configure();

    // Get the Xios singleton instance and check calendar step is zero initially
    // NOTE: The singleton is created when Xios::getInstance() is first called. In this test, this
    //       happens when the time sets set by ModelMetadata::setTime(). This occurs in the call to
    //       Model::configure() above.
    Xios& xiosHandler = Xios::getInstance();
    REQUIRE(xiosHandler.getCalendarStep() == 0);

    // Set ModelArray dimensions
    REQUIRE(ModelArray::nComponents(ModelArray::Type::DG) == DGCOMP);
    REQUIRE(ModelArray::size(ModelArray::Dimension::DG) == DGCOMP);

    // Read initial state from the test input file, checking that the it contains the expected
    // fields
    // NOTE: XIOS axes, domains, and grids are created by the ParaGridIO constructor, which is
    //       constructed in the call to StructureFactory::stateFromFile(). The XIOS context also
    //       gets closed in this call.
    ModelState modelstate = StructureFactory::stateFromFile(inputFilename);
    REQUIRE(modelstate.data.count(maskName) > 0);
    REQUIRE(modelstate.data.count(longitudeName) > 0);
    REQUIRE(modelstate.data.count(latitudeName) > 0);
    REQUIRE(modelstate.data.count(gridAzimuthName) > 0);
    REQUIRE(modelstate.data.count(coordsName) > 0);
    REQUIRE(modelstate.data.count(ciceName) > 0);
    REQUIRE(modelstate.data.count(hiceName) > 0);
    REQUIRE(modelstate.data.count(hsnowName) > 0);
    REQUIRE(modelstate.data.count(ticeName) > 0);
    REQUIRE(modelstate.data.count(shearName) > 0);

    // Check files with the expected names don't exist yet
    REQUIRE_FALSE(std::filesystem::exists("readwrite*.nc"));

    // Simulate 4 iterations (timesteps)
    ModelMetadata& metadata = ModelMetadata::getInstance();
    const Duration& timestep = metadata.stepLength();
    int rank;
    MPI_Comm_rank(test_comm, &rank);
    for (int ts = 0; ts <= 4; ts++) {
        StructureFactory::fileFromState(modelstate, outputFilename, true);

        // Update the current timestep and verify it's updated in XIOS
        metadata.incrementTime(timestep);
        REQUIRE(xiosHandler.getCalendarStep() == ts + 1);
    }

    // Check the files have indeed been created then remove them
    REQUIRE(std::filesystem::exists("readwrite_2023-03-17T17:11:00Z-2023-03-17T20:10:59Z.nc"));
    REQUIRE(std::filesystem::exists("readwrite_2023-03-17T20:11:00Z-2023-03-17T23:10:59Z.nc"));
    if (rank == 0) {
        std::filesystem::remove("readwrite_2023-03-17T17:11:00Z-2023-03-17T20:10:59Z.nc");
        std::filesystem::remove("readwrite_2023-03-17T20:11:00Z-2023-03-17T23:10:59Z.nc");
    }

    xiosHandler.context_finalize();
    Finalizer::finalize();
}

}
