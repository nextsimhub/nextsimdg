/*!
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @brief   Tests for XIOS functionality for writing restart files.
 * @details The functionality of writing restarts via XIOS is tested.
 */
#include <doctest/extensions/doctest_mpi.h>
#undef INFO

#include "StructureModule/include/ParametricGrid.hpp"
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
const std::string outputFilename = testFilesDir + "/restart%Y-%m-%dT%H:%M:%SZ.nc";

static const int DGCOMP = 6;
static const int DGSTRESSCOMP = 8;
static const int CGDEGREE = 2;

namespace Nextsim {

/*!
 * TestXiosWriteRestart
 *
 * Test writing of restarts via `dumpModelState`.
 */
MPI_TEST_CASE("TestXiosWriteRestart", 2)
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
    config << "[XiosOutput]" << std::endl;
    config << "field_names = " << maskName << "," << longitudeName << "," << latitudeName << ","
           << gridAzimuthName << "," << ciceName << "," << hiceName << "," << damageName << ","
           << hsnowName << "," << ticeName << "," << uName << "," << std::endl;
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

    // Set ModelArray dimensions
    const size_t nx = ModelArray::size(ModelArray::Dimension::X);
    const size_t ny = ModelArray::size(ModelArray::Dimension::Y);
    REQUIRE(ModelArray::nComponents(ModelArray::Type::DG) == DGCOMP);
    REQUIRE(ModelArray::size(ModelArray::Dimension::DG) == DGCOMP);

    // The ParametricGrid structure is required by XIOS
    Module::setImplementation<IStructure>("Nextsim::ParametricGrid");

    // Set field types for restarts
    xiosHandler.setFieldType(longitudeName, ModelArray::Type::H);
    xiosHandler.setFieldType(latitudeName, ModelArray::Type::H);
    xiosHandler.setFieldType(maskName, ModelArray::Type::H);
    xiosHandler.setFieldType(ciceName, ModelArray::Type::DG);
    xiosHandler.setFieldType(hiceName, ModelArray::Type::DG);
    xiosHandler.setFieldType(damageName, ModelArray::Type::DG);
    xiosHandler.setFieldType(hsnowName, ModelArray::Type::DG);
    xiosHandler.setFieldType(gridAzimuthName, ModelArray::Type::VERTEX);
    xiosHandler.setFieldType(ticeName, ModelArray::Type::DGSTRESS);
    xiosHandler.setFieldType(uName, ModelArray::Type::CG);

    // Create some fake data to test writing methods
    HField longitude(ModelArray::Type::H);
    longitude.resize();
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            longitude(i, j) = i;
        }
    }
    HField latitude(ModelArray::Type::H);
    latitude.resize();
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            latitude(i, j) = j;
        }
    }
    /*
     * Mask definition, where 0 indicates land and 1 indicates ocean:
     *
     * Rank 0:  Rank 1:
     * -----    -----
     * |0|1|    |0|1|
     * -----    -----
     * |1|1|    |1|1|
     * -----    -----
     *
     * That is, mask is zero when i = 0 and j = 0 and one otherwise.
     */
    HField mask(ModelArray::Type::H);
    mask.resize();
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            mask(i, j) = (i == 0 && j == 0 ? 0.0 : 1.0);
        }
    }
    DGField cice(ModelArray::Type::DG);
    cice.resize();
    DGField hice(ModelArray::Type::DG);
    hice.resize();
    DGField damage(ModelArray::Type::DG);
    damage.resize();
    DGField hsnow(ModelArray::Type::DG);
    hsnow.resize();
    VertexField grid_azimuth(ModelArray::Type::VERTEX);
    grid_azimuth.resize();
    DGSField tice(ModelArray::Type::DGSTRESS);
    tice.resize();
    CGField uice(ModelArray::Type::CG);
    uice.resize();

    // Check files with the expected names don't exist yet
    REQUIRE_FALSE(std::filesystem::exists("restart*.nc"));

    // Check calendar step is zero initially
    REQUIRE(xiosHandler.getCalendarStep() == 0);

    // Simulate 4 iterations (timesteps)
    ModelMetadata& metadata = ModelMetadata::getInstance();
    const Duration& timestep = metadata.stepLength();
    int rank;
    MPI_Comm_rank(test_comm, &rank);
    for (int ts = 0; ts <= 4; ts++) {

        // Update DGField restarts
        // NOTE: NaN values for mask when i = 0 and j = 0
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {
                for (size_t d = 0; d < DGCOMP; ++d) {
                    const float value
                        = (i == 0 && j == 0) ? NAN : 1.0 * ts * (d + DGCOMP * (i + nx * j));
                    cice.components({ i, j })[d] = value;
                    hice.components({ i, j })[d] = value;
                    damage.components({ i, j })[d] = value;
                    hsnow.components({ i, j })[d] = value;
                }
            }
        }

        // Update VertexField restarts
        for (size_t j = 0; j < ny + 1; ++j) {
            for (size_t i = 0; i < nx + 1; ++i) {
                if (rank == 0) {
                    grid_azimuth.components({ i, j })[0] = 1.0 * ts * i;
                    grid_azimuth.components({ i, j })[1] = 1.0 * ts * j;
                } else {
                    grid_azimuth.components({ i, j })[0] = 1.0 * ts * (i + 2);
                    grid_azimuth.components({ i, j })[1] = 1.0 * ts * j;
                }
            }
        }

        // Update DGSField restarts
        // NOTE: NaN values for mask when i = 0 and j = 0
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {
                for (size_t d = 0; d < DGSTRESSCOMP; ++d) {
                    const float value
                        = (i == 0 && j == 0) ? NAN : 2.0 * ts * (d + DGSTRESSCOMP * (i + nx * j));
                    tice.components({ i, j })[d] = value;
                }
            }
        }

        /*
         * Update CGField restarts
         *
         * In this test, the CG field has degree 2, so there are 9 DoFs associated with each cell,
         * although many of these are shared with neighbouring cells.
         *
         * On rank 0, the top-left cell is land, so the nodal values are masked except for those
         * shared with ocean cells, which effectively become boundary nodes.
         *
         * On rank 1, we also need to account for the the shared nodes on the left-hand-side between
         * the two subdomains, which must take consistent values on either side. Again, these
         * effectively become boundary nodes.
         *
         * Rank 0:    Rank 1:
         * 0-0-1-1-1  1-0-1-1-1
         * | | | | |  | | | | |
         * 0-0-1-1-1  1-0-1-1-1
         * | | | | |  | | | | |
         * 1-1-1-1-1  1-1-1-1-1
         * | | | | |  | | | | |
         * 1-1-1-1-1  1-1-1-1-1
         * | | | | |  | | | | |
         * 1-1-1-1-1  1-1-1-1-1
         *
         * That is, mask is zero when i <= 1 and j <= 1 on rank 0 and when i = j = 1 on rank 1.
         */
        for (size_t j = 0; j < CGDEGREE * ny + 1; ++j) {
            for (size_t i = 0; i < CGDEGREE * nx + 1; ++i) {
                float value;
                if (rank == 0) {
                    value = (i <= 1 && j <= 1) ? NAN : 1.0 * ts * ((i + 1) * (j + 1));
                } else {
                    value = (i == 1 && j == 1) ? NAN : 1.0 * ts * ((i + 5) * (j + 1));
                }
                uice(i, j) = value;
            }
        }

        // Set up ModelState for restarts and write out
        ModelState restarts = { {
                                    { maskName, mask },
                                    { longitudeName, longitude },
                                    { latitudeName, latitude },
                                    { ciceName, cice },
                                    { hiceName, hice },
                                    { damageName, damage },
                                    { hsnowName, hsnow },
                                    { gridAzimuthName, grid_azimuth },
                                    { ticeName, tice },
                                    { uName, uice },
                                },
            {} };
        StructureFactory::fileFromState(restarts, outputFilename, true);

        // Update the current timestep and verify it's updated in XIOS
        metadata.incrementTime(timestep);
        REQUIRE(xiosHandler.getCalendarStep() == ts + 1);
    }

    // Check the files have indeed been created
    // NOTE: Don't remove them because their contents are checked in XiosReadRestart_test
    REQUIRE(std::filesystem::exists("restart_2023-03-17T17:11:00Z-2023-03-17T20:10:59Z.nc"));
    REQUIRE(std::filesystem::exists("restart_2023-03-17T20:11:00Z-2023-03-17T23:10:59Z.nc"));

    xiosHandler.context_finalize();
    Finalizer::finalize();
}

}
