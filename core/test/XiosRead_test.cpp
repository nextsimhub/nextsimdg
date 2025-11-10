/*!
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @brief   Tests for XIOS read functionality
 * @details
 * This test is designed to test the file reading functionality of the C++
 * interface for XIOS.
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
const std::string restartFilename = testFilesDir + "/xios_test_input.nc";
const std::string forcingFilename = testFilesDir + "/xios_test_forcing.nc";

static const int DGCOMP = 6;
static const int DGSTRESSCOMP = 8;
static const int CGDEGREE = 2;

namespace Nextsim {

/*!
 * TestXiosRead
 *
 * This function tests the file reading functionality of the C++ interface for XIOS for fields with
 * two and three spatial dimensions. The test runs with two MPI ranks.
 */
MPI_TEST_CASE("TestXiosRead", 2)
{
    std::stringstream config;
    config << "[model]" << std::endl;
    config << "start = 2023-03-17T17:11:00Z" << std::endl;
    config << "stop = 2023-03-17T23:11:00Z" << std::endl;
    config << "time_step = P0-0T01:30:00" << std::endl;
    config << "init_file = " << restartFilename << std::endl;
    config << "restart_period = P0-0T01:30:00" << std::endl;
    config << "partition_file = xios_test_partition_metadata_2.nc" << std::endl;
    config << "[XiosInput]" << std::endl;
    config << "field_names = " << maskName << "," << coordsName << "," << hiceName << ","
           << ticeName << "," << uName << std::endl;
    config << "[XiosForcing]" << std::endl;
    config << "filename = " << forcingFilename << std::endl;
    config << "field_names = " << hsnowName << std::endl;
    config << "period = P0-0T01:30:00" << std::endl;
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

    // Check the input files exists
    REQUIRE(std::filesystem::exists(restartFilename));
    REQUIRE(std::filesystem::exists(forcingFilename));

    // Check calendar step is zero initially
    REQUIRE(xiosHandler.getCalendarStep() == 0);

    // Deduce the local lengths of the two dimensions
    const size_t nx = ModelArray::size(ModelArray::Dimension::X);
    const size_t ny = ModelArray::size(ModelArray::Dimension::Y);
    REQUIRE(ModelArray::nComponents(ModelArray::Type::DG) == DGCOMP);
    REQUIRE(ModelArray::size(ModelArray::Dimension::DG) == DGCOMP);

    // Read restarts from file and check they take the expected values
    ModelState restarts = grid.getModelState(restartFilename);
    int rank;
    MPI_Comm_rank(test_comm, &rank);
    for (auto& entry : restarts.data) {
        if (entry.first == maskName) {
            for (size_t j = 0; j < ny; ++j) {
                for (size_t i = 0; i < nx; ++i) {
                    REQUIRE(entry.second(i, j) == doctest::Approx(j >= 1 ? 1.0 : 0.0));
                }
            }
        } else if (entry.first == coordsName) {
            for (size_t j = 0; j < ny + 1; ++j) {
                for (size_t i = 0; i < nx + 1; ++i) {
                    if (rank == 0) {
                        REQUIRE(entry.second.components({ i, j })[0] == doctest::Approx(i));
                    } else {
                        REQUIRE(entry.second.components({ i, j })[0] == doctest::Approx(i + 2));
                    }
                    REQUIRE(entry.second.components({ i, j })[1] == doctest::Approx(j));
                }
            }
        } else if (entry.first == hiceName) {
            for (size_t j = 0; j < ny; ++j) {
                for (size_t i = 0; i < nx; ++i) {
                    for (size_t d = 0; d < DGCOMP; ++d) {
                        float expected = 1.0 * (d + DGCOMP * (i + nx * j));
                        REQUIRE(entry.second.components({ i, j })[d] == doctest::Approx(expected));
                    }
                }
            }
        } else if (entry.first == ticeName) {
            for (size_t j = 0; j < ny; ++j) {
                for (size_t i = 0; i < nx; ++i) {
                    for (size_t d = 0; d < DGSTRESSCOMP; ++d) {
                        float expected = 2.0 * (d + DGSTRESSCOMP * (i + nx * j));
                        REQUIRE(entry.second.components({ i, j })[d] == doctest::Approx(expected));
                    }
                }
            }
        } else if (entry.first == uName) {
            for (size_t j = 0; j < CGDEGREE * ny + 1; ++j) {
                for (size_t i = 0; i < CGDEGREE * nx + 1; ++i) {
                    if (rank == 0) {
                        REQUIRE(entry.second(i, j) == doctest::Approx((i + 1) * (j + 1)));
                    } else {
                        REQUIRE(entry.second(i, j) == doctest::Approx((i + 5) * (j + 1)));
                    }
                }
            }
        }
    }

    // Simulate 4 iterations (timesteps), reading forcing data at each
    ModelMetadata& metadata = ModelMetadata::getInstance();
    Duration timestep = metadata.stepLength();
    // TODO: Avoid making configGetForcingFieldNames public?
    auto forcingFieldNames = xiosHandler.configGetForcingFieldNames();
    for (int ts = 0; ts <= 4; ts++) {

        // Read forcings from file and check they take the expected values
        TimePoint time = xiosHandler.getCurrentDate();
        ModelState forcings = pio->readForcingTimeStatic(forcingFieldNames, time, forcingFilename);
        for (auto& entry : forcings.data) {
            REQUIRE(entry.first == hsnowName);
            for (size_t j = 0; j < ny; ++j) {
                for (size_t i = 0; i < nx; ++i) {
                    REQUIRE(entry.second(i, j) == doctest::Approx(0.1 * ts));
                }
            }
        }

        // Update the current timestep and verify it's updated in XIOS
        metadata.incrementTime(timestep);
        REQUIRE(xiosHandler.getCalendarStep() == ts + 1);
    }

    xiosHandler.context_finalize();
    Finalizer::finalize();
}
}
