/*!
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @brief   Tests for XIOS functionality for reading restart files.
 * @details The functionality of reading restarts via XIOS is tested.
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
const std::string restartFilename
    = testFilesDir + "/restart_2023-03-17T17:11:00Z-2023-03-17T20:10:59Z.nc";

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
    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    // Check the input file exists
    if (!std::filesystem::exists(restartFilename)) {
        throw std::runtime_error(
            "XiosReadRestart_test: Input file not found. Did you run XiosWrite_test?");
    }

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

    // Deduce the local lengths of the two dimensions
    const size_t nx = ModelArray::size(ModelArray::Dimension::X);
    const size_t ny = ModelArray::size(ModelArray::Dimension::Y);
    REQUIRE(ModelArray::nComponents(ModelArray::Type::DG) == DGCOMP);
    REQUIRE(ModelArray::size(ModelArray::Dimension::DG) == DGCOMP);

    // Read restarts from file and check they take the expected values
    // NOTE: The ParametricGrid is created and the XIOS context definition is closed in the call to
    //       StructureFactory::stateFromFile()
    int rank;
    MPI_Comm_rank(test_comm, &rank);
    float ts = 2; // Corresponds to 2023-03-17T20:10:59Z
    ModelState modelstate = StructureFactory::stateFromFile(restartFilename);
    ModelArray& mask = modelstate.data.at(maskName);
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            REQUIRE(mask(i, j) == doctest::Approx(i == 0 && j == 0 ? 0.0 : 1.0));
        }
    }
    for (const auto& [fieldName, modelarray] : modelstate.data) {
        if (fieldName == longitudeName) {
            for (size_t j = 0; j < ny; ++j) {
                for (size_t i = 0; i < nx; ++i) {
                    REQUIRE(modelarray(i, j) == doctest::Approx(i));
                }
            }
        } else if (fieldName == latitudeName) {
            for (size_t j = 0; j < ny; ++j) {
                for (size_t i = 0; i < nx; ++i) {
                    REQUIRE(modelarray(i, j) == doctest::Approx(j));
                }
            }
        } else if (fieldName == gridAzimuthName) {
            for (size_t j = 0; j < ny; ++j) {
                for (size_t i = 0; i < nx; ++i) {
                    REQUIRE(modelarray(i, j) == doctest::Approx(0.0));
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
        } else if (fieldName == ciceName || fieldName == hiceName || fieldName == damageName
            || fieldName == hsnowName) {
            for (size_t j = 0; j < ny; ++j) {
                for (size_t i = 0; i < nx; ++i) {
                    for (size_t d = 0; d < DGCOMP; ++d) {
                        if (i == 0 && j == 0) {
                            REQUIRE(std::isnan(modelarray.components({ i, j })[d]));
                        } else {
                            const float expected = ts * (d + DGCOMP * (i + nx * j));
                            REQUIRE(
                                modelarray.components({ i, j })[d] == doctest::Approx(expected));
                        }
                    }
                }
            }
        } else if (fieldName == ticeName) {
            for (size_t j = 0; j < ny; ++j) {
                for (size_t i = 0; i < nx; ++i) {
                    for (size_t d = 0; d < DGSTRESSCOMP; ++d) {
                        if (i == 0 && j == 0) {
                            REQUIRE(std::isnan(modelarray.components({ i, j })[d]));
                        } else {
                            const float expected = 2.0 * ts * (d + DGSTRESSCOMP * (i + nx * j));
                            REQUIRE(
                                modelarray.components({ i, j })[d] == doctest::Approx(expected));
                        }
                    }
                }
            }
        } else if (fieldName == shearName) {
            for (size_t j = 0; j < CGDEGREE * ny + 1; ++j) {
                for (size_t i = 0; i < CGDEGREE * nx + 1; ++i) {
                    if (rank == 0) {
                        if (i <= 1 && j <= 1) {
                            REQUIRE(std::isnan(modelarray(i, j)));
                        } else {
                            REQUIRE(modelarray(i, j) == doctest::Approx(ts * (i + 1) * (j + 1)));
                        }
                    } else {
                        if (i == 1 && j == 1) {
                            REQUIRE(std::isnan(modelarray(i, j)));
                        } else {
                            REQUIRE(modelarray(i, j) == doctest::Approx(ts * (i + 5) * (j + 1)));
                        }
                    }
                }
            }
        }
    }

    // Remove the restart files
    if (rank == 0) {
        std::filesystem::remove("restart_2023-03-17T17:11:00Z-2023-03-17T20:10:59Z.nc");
        std::filesystem::remove("restart_2023-03-17T20:11:00Z-2023-03-17T23:10:59Z.nc");
    }

    xiosHandler.context_finalize();
    Finalizer::finalize();
}
}
