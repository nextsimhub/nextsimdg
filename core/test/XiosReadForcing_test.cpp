/*!
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @brief   Tests for XIOS functionality for reading forcings
 * @details The functionality of reading forcings via XIOS is tested.
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
const std::string inputFilename = testFilesDir + "/xios_test_input.nc";
const std::string forcingFilename = testFilesDir + "/xios_test_forcing.nc";

static const int DGCOMP = 6;
static const int DGSTRESSCOMP = 8;
static const int CGDEGREE = 2;

namespace Nextsim {

/*!
 * TestXiosReadForcing
 *
 * Test reading of restarts via `readForcingTimeStatic`.
 */
MPI_TEST_CASE("TestXiosReadForcing", 2)
{
    std::stringstream config;
    config << "[model]" << std::endl;
    config << "start = 2023-03-17T17:11:00Z" << std::endl;
    config << "stop = 2023-03-17T23:11:00Z" << std::endl;
    config << "time_step = P0-0T01:30:00" << std::endl;
    config << "init_file = " << inputFilename << std::endl;
    config << "restart_period = P0-0T06:00:00" << std::endl;
    config << "partition_file = xios_test_partition_metadata_2.nc" << std::endl;
    config << "[XiosForcing]" << std::endl;
    config << "filename = " << forcingFilename << std::endl;
    config << "field_names = " << hsnowName << std::endl;
    config << "period = P0-0T03:00:00" << std::endl;
    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    // Create ModelMPI instance based off the test communicator
    auto& modelMPI = ModelMPI::getInstance(test_comm);

    // Create a Model and configure it so that time options are parsed
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

    // NOTE: Needs calling before Xios::getCurrentDate()
    xiosHandler.close_context_definition();

    // Check the input file exists
    REQUIRE(std::filesystem::exists(forcingFilename));

    // Check calendar step is zero initially
    REQUIRE(xiosHandler.getCalendarStep() == 0);

    // Deduce the local lengths of the two dimensions
    const size_t nx = ModelArray::size(ModelArray::Dimension::X);
    const size_t ny = ModelArray::size(ModelArray::Dimension::Y);
    REQUIRE(ModelArray::nComponents(ModelArray::Type::DG) == DGCOMP);
    REQUIRE(ModelArray::size(ModelArray::Dimension::DG) == DGCOMP);

    // Simulate 4 iterations (timesteps), reading forcing data at each
    ModelMetadata& metadata = ModelMetadata::getInstance();
    const Duration& timestep = metadata.stepLength();
    for (int ts = 0; ts <= 4; ts += 2) {

        // Read forcings from file and check they take the expected values
        // TODO: Avoid making forcingFieldNames public?
        const TimePoint& time = xiosHandler.getCurrentDate();
        const ModelState forcings
            = pio->readForcingTimeStatic(xiosHandler.forcingFieldNames, time, forcingFilename);
        for (const auto& [fieldName, modelarray] : forcings.data) {
            REQUIRE(fieldName == hsnowName);
            for (size_t j = 0; j < ny; ++j) {
                for (size_t i = 0; i < nx; ++i) {
                    REQUIRE(modelarray(i, j) == doctest::Approx(0.1 * ts));
                }
            }
        }

        // Update the current timestep and verify it's updated in XIOS
        metadata.incrementTime(timestep);
        metadata.incrementTime(timestep);
        REQUIRE(xiosHandler.getCalendarStep() == ts + 2);
    }

    xiosHandler.context_finalize();
    Finalizer::finalize();
}
}
