/*!
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @brief   Tests for XIOS functionality for reading restart files.
 * @details The functionality of reading restarts via XIOS is tested.
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
const std::string restartFilename
    = testFilesDir + "/xios_test_restart_20230317171100-20230317201059.nc";

static const int DGCOMP = 6;
static const int DGSTRESSCOMP = 8;
static const int CGDEGREE = 2;

namespace Nextsim {

/*!
 * TestXiosReadRestart
 *
 * Test reading of restarts via `getModelState`.
 */
MPI_TEST_CASE("TestXiosReadRestart", 2)
{
    std::stringstream config;
    config << "[model]" << std::endl;
    config << "start = 2023-03-17T17:11:00Z" << std::endl;
    config << "stop = 2023-03-17T23:11:00Z" << std::endl;
    config << "time_step = P0-0T01:30:00" << std::endl;
    config << "init_file = " << restartFilename << std::endl;
    config << "restart_period = P0-0T03:00:00" << std::endl;
    config << "partition_file = xios_test_partition_metadata_2.nc" << std::endl;
    config << "[XiosInput]" << std::endl;
    config << "field_names = " << maskName << "," << coordsName << "," << hiceName << ","
           << ticeName << "," << uName << std::endl;
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
    if (!std::filesystem::exists(restartFilename)) {
        throw std::runtime_error(
            "XiosReadRestart_test: Input file not found. Did you run XiosWrite_test?");
    }

    // Check calendar step is zero initially
    REQUIRE(xiosHandler.getCalendarStep() == 0);

    // Deduce the local lengths of the two dimensions
    const size_t nx = ModelArray::size(ModelArray::Dimension::X);
    const size_t ny = ModelArray::size(ModelArray::Dimension::Y);
    REQUIRE(ModelArray::nComponents(ModelArray::Type::DG) == DGCOMP);
    REQUIRE(ModelArray::size(ModelArray::Dimension::DG) == DGCOMP);

    // Read restarts from file and check they take the expected values
    int rank;
    MPI_Comm_rank(test_comm, &rank);
    float ts = 2; // Corresponds to 2023-03-17T20:10:59Z
    for (const auto [fieldName, modelarray] : grid.getModelState(restartFilename).data) {
        if (fieldName == maskName) {
            for (size_t j = 0; j < ny; ++j) {
                for (size_t i = 0; i < nx; ++i) {
                    REQUIRE(modelarray(i, j) == doctest::Approx(j >= 1 ? 1.0 : 0.0));
                }
            }
        } else if (fieldName == coordsName) {
            for (size_t j = 0; j < ny + 1; ++j) {
                for (size_t i = 0; i < nx + 1; ++i) {
                    float expected_x;
                    if (rank == 0) {
                        expected_x = ts * i;
                    } else {
                        expected_x = ts * (i + 2);
                    }
                    const float expected_y = ts * j;
                    REQUIRE(modelarray.components({ i, j })[0] == doctest::Approx(expected_x));
                    REQUIRE(modelarray.components({ i, j })[1] == doctest::Approx(expected_y));
                }
            }
        } else if (fieldName == hiceName) {
            for (size_t j = 0; j < ny; ++j) {
                for (size_t i = 0; i < nx; ++i) {
                    for (size_t d = 0; d < DGCOMP; ++d) {
                        const float expected = ts * (d + DGCOMP * (i + nx * j));
                        REQUIRE(modelarray.components({ i, j })[d] == doctest::Approx(expected));
                    }
                }
            }
        } else if (fieldName == ticeName) {
            for (size_t j = 0; j < ny; ++j) {
                for (size_t i = 0; i < nx; ++i) {
                    for (size_t d = 0; d < DGSTRESSCOMP; ++d) {
                        const float expected = 2.0 * ts * (d + DGSTRESSCOMP * (i + nx * j));
                        REQUIRE(modelarray.components({ i, j })[d] == doctest::Approx(expected));
                    }
                }
            }
        } else if (fieldName == uName) {
            for (size_t j = 0; j < CGDEGREE * ny + 1; ++j) {
                for (size_t i = 0; i < CGDEGREE * nx + 1; ++i) {
                    float expected;
                    if (rank == 0) {
                        expected = ts * (i + 1) * (j + 1);
                    } else {
                        expected = ts * (i + 5) * (j + 1);
                    }
                    REQUIRE(modelarray(i, j) == doctest::Approx(expected));
                }
            }
        }
    }

    // Remove the restart files
    if (rank == 0) {
        std::filesystem::remove("xios_test_restart_20230317171100-20230317201059.nc");
        std::filesystem::remove("xios_test_restart_20230317201100-20230317231059.nc");
    }

    xiosHandler.context_finalize();
    Finalizer::finalize();
}
}
