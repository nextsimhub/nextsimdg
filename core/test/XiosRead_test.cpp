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

static const int DG = 3;
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
    // Enable XIOS in the 'config' and provide parameters to configure it
    enableXios();
    std::stringstream config;
    config << "[model]" << std::endl;
    config << "start = 2023-03-17T17:11:00Z" << std::endl;
    config << "stop = 2023-03-17T23:11:00Z" << std::endl;
    config << "time_step = P0-0T01:30:00" << std::endl;
    config << "init_file = " << restartFilename << std::endl;
    config << "restart_period = P0-0T01:30:00" << std::endl;
    config << "[XiosInput]" << std::endl;
    config << "field_names = " << maskName << "," << coordsName << "," << hiceName << ","
           << ticeName << "," << uName << std::endl;
    config << "[XiosForcing]" << std::endl;
    config << "filename = " << forcingFilename << std::endl;
    config << "field_names = " << hsnowName << std::endl;
    config << "period = P0-0T01:30:00" << std::endl;
    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    // Create ModelMetadata instance based off a partition metadata file
    auto& modelMPI = ModelMPI::getInstance(test_comm);
    auto& metadata = ModelMetadata::getInstance("xios_test_partition_metadata_2.nc");

    // Create a Model and configure it so that time options are parsed
    Model model;
    model.configureTime(); // TODO: Use Model.configure to parse restart files this way, too?

    // Create ParametricGrid and ParaGridIO instances
    Module::setImplementation<IStructure>("Nextsim::ParametricGrid");
    ParametricGrid grid;
    ParaGridIO* pio = new ParaGridIO(grid);
    grid.setIO(pio);

    // TODO: We could deduce this from the NetCDF file
    ModelArray::setNComponents(ModelArray::Type::DG, DG);
    ModelArray::setNComponents(ModelArray::Type::DGSTRESS, DGSTRESSCOMP);
    ModelArray::setNComponents(ModelArray::Type::VERTEX, ModelArray::nCoords);
    REQUIRE(ModelArray::nComponents(ModelArray::Type::DG) == DG);
    REQUIRE(ModelArray::nComponents(ModelArray::Type::DGSTRESS) == DGSTRESSCOMP);
    REQUIRE(ModelArray::nComponents(ModelArray::Type::VERTEX) == ModelArray::nCoords);

    // Get the Xios singleton instance and check it's initialized
    Xios& xiosHandler = Xios::getInstance();
    REQUIRE(xiosHandler.isInitialized());

    // Affix ModelMetadata to Xios handler
    // TODO: Automate this - can't be inlined in Xios::getInstance because need set field types
    xiosHandler.affixModelMetadata();

    xiosHandler.close_context_definition();

    // Check calendar step is zero initially
    REQUIRE(xiosHandler.getCalendarStep() == 0);

    // Check the input files exists
    REQUIRE(std::filesystem::exists(restartFilename));
    REQUIRE(std::filesystem::exists(forcingFilename));

    // Check calendar step is zero initially
    metadata.setTime(xiosHandler.getCalendarStart());
    REQUIRE(xiosHandler.getCalendarStep() == 0);

    // Deduce the local lengths of the two dimensions
    const size_t nx = ModelArray::size(ModelArray::Dimension::X);
    const size_t ny = ModelArray::size(ModelArray::Dimension::Y);

    // Read restarts from file and check they take the expected values
    metadata.setTime(xiosHandler.getCalendarStart());
    REQUIRE(xiosHandler.getCalendarStep() == 0);
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
                    for (size_t d = 0; d < DG; ++d) {
                        float expected = 1.0 * (d + DG * (i + nx * j));
                        REQUIRE(entry.second.components({ i, j })[d] == doctest::Approx(expected));
                    }
                }
            }
        } else if (entry.first == ticeName) {
            for (size_t j = 0; j < ny; ++j) {
                for (size_t i = 0; i < nx; ++i) {
                    for (size_t d = 0; d < DGSTRESSCOMP; ++d) {
                        REQUIRE(entry.second.components({ i, j })[d]
                            == doctest::Approx(2.0 * (d + DGSTRESSCOMP * (i + nx * j))));
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
    Duration timestep = xiosHandler.getCalendarTimestep();
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
