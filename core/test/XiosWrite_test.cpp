/*!
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @brief   Tests for XIOS write functionality
 * @details
 * This test is designed to test the file writing functionality of the C++
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
const std::string restartFilename = testFilesDir + "/xios_test_output.nc";
const std::string diagnosticFilename = testFilesDir + "/xios_test_diagnostic.nc";

static const int DG = 3;
static const int DGSTRESSCOMP = 8;
static const int CGDEGREE = 2;

namespace Nextsim {

/*!
 * TestXiosWrite
 *
 * This function tests the file writing functionality of the C++ interface for XIOS for fields with
 * two and three spatial dimensions. The test runs with two MPI ranks.
 */
MPI_TEST_CASE("TestXiosWrite", 2)
{
    // Enable XIOS in the 'config' and provide parameters to configure it
    enableXios();
    std::stringstream config;
    config << "[model]" << std::endl;
    config << "start = 2023-03-17T17:11:00Z" << std::endl;
    config << "stop = 2023-03-17T23:11:00Z" << std::endl;
    config << "time_step = P0-0T01:30:00" << std::endl;
    config << "restart_file = " << restartFilename << std::endl;
    config << "restart_period = P0-0T01:30:00" << std::endl;
    config << "[XiosOutput]" << std::endl;
    config << "field_names = " << maskName << "," << coordsName << "," << hiceName << ","
           << ticeName << "," << ciceName << std::endl;
    config << "period = P0-0T01:30:00" << std::endl;
    config << "[XiosDiagnostic]" << std::endl;
    config << "filename = " << diagnosticFilename << std::endl;
    config << "field_names = " << uName << std::endl;
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

    // Get the Xios singleton instance and check it's initialized
    Xios& xiosHandler = Xios::getInstance();
    REQUIRE(xiosHandler.isInitialized());
    REQUIRE(xiosHandler.getClientMPISize() == 2);

    // Set ModelArray dimensions
    const size_t nx_glo = 4;
    const size_t ny_glo = 2;
    const size_t nx = 2;
    const size_t ny = 2;
    size_t nx_start;
    const size_t ny_start = 0;
    const int rank = xiosHandler.getClientMPIRank();
    if (rank == 0) {
        nx_start = 0;
    } else {
        nx_start = nx_glo - nx;
    }
    ModelArray::setDimension(ModelArray::Dimension::X, nx_glo, nx, nx_start);
    ModelArray::setDimension(ModelArray::Dimension::XVERTEX, nx_glo + 1, nx + 1, nx_start);
    ModelArray::setDimension(
        ModelArray::Dimension::XCG, CGDEGREE * nx_glo + 1, CGDEGREE * nx + 1, CGDEGREE * nx_start);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny_glo, ny, ny_start);
    ModelArray::setDimension(ModelArray::Dimension::YVERTEX, ny_glo + 1, ny + 1, ny_start);
    ModelArray::setDimension(
        ModelArray::Dimension::YCG, CGDEGREE * ny_glo + 1, CGDEGREE * ny + 1, ny_start);
    ModelArray::setNComponents(ModelArray::Type::DG, DG);
    ModelArray::setNComponents(ModelArray::Type::DGSTRESS, DGSTRESSCOMP);
    ModelArray::setNComponents(ModelArray::Type::VERTEX, ModelArray::nCoords);
    REQUIRE(ModelArray::nComponents(ModelArray::Type::DG) == DG);
    REQUIRE(ModelArray::nComponents(ModelArray::Type::DGSTRESS) == DGSTRESSCOMP);
    REQUIRE(ModelArray::nComponents(ModelArray::Type::VERTEX) == ModelArray::nCoords);

    // Affix ModelMetadata to Xios handler
    // TODO: Automate this - can't be inlined in Xios::getInstance because need set field types
    xiosHandler.affixModelMetadata();

    // Set field types
    xiosHandler.setFieldType(maskName, ModelArray::Type::H);
    xiosHandler.setFieldType(coordsName, ModelArray::Type::VERTEX);
    xiosHandler.setFieldType(hiceName, ModelArray::Type::DG);
    xiosHandler.setFieldType(ticeName, ModelArray::Type::DGSTRESS);
    xiosHandler.setFieldType(ciceName, ModelArray::Type::CG);
    xiosHandler.setFieldType(uName, ModelArray::Type::H);

    // Set file split frequency for restarts (but not diagnostics)
    // NOTE: Files are created when the XIOS handler is constructed
    const std::string fileId = "xios_test_output";
    xiosHandler.setFileSplitFreq(fileId, Duration("P0-0T03:00:00"));

    xiosHandler.close_context_definition();

    // Create some fake data to test writing methods
    HField mask(ModelArray::Type::H);
    mask.resize();
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            mask(i, j) = j >= 1 ? 1.0 : 0.0;
        }
    }
    VertexField coordinates(ModelArray::Type::VERTEX);
    coordinates.resize();
    for (size_t j = 0; j < ny + 1; ++j) {
        for (size_t i = 0; i < nx + 1; ++i) {
            if (rank == 0) {
                coordinates.components({ i, j })[0] = (double)i;
                coordinates.components({ i, j })[1] = (double)j;
            } else {
                coordinates.components({ i, j })[0] = (double)(i + 2);
                coordinates.components({ i, j })[1] = (double)j;
            }
        }
    }
    DGField hice(ModelArray::Type::DG);
    hice.resize();
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            for (size_t d = 0; d < DG; ++d) {
                hice.components({ i, j })[d] = 1.0 * (d + DG * (i + nx * j));
            }
        }
    }
    DGSField tice(ModelArray::Type::DGSTRESS);
    tice.resize();
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            for (size_t d = 0; d < DGSTRESSCOMP; ++d) {
                tice.components({ i, j })[d] = 2.0 * (d + DGSTRESSCOMP * (i + nx * j));
            }
        }
    }
    CGField cice(ModelArray::Type::CG);
    cice.resize();
    for (size_t j = 0; j < CGDEGREE * ny + 1; ++j) {
        for (size_t i = 0; i < CGDEGREE * nx + 1; ++i) {
            if (rank == 0) {
                cice(i, j) = (double)((i + 1) * (j + 1));
            } else {
                cice(i, j) = (double)((i + 5) * (j + 1));
            }
        }
    }
    HField u(ModelArray::Type::H);
    u.resize();

    // Check files with the expected names don't exist yet
    REQUIRE_FALSE(std::filesystem::exists("xios_test_output*.nc"));
    REQUIRE_FALSE(std::filesystem::exists("xios_test_diagnostic*.nc"));

    // Simulate 4 iterations (timesteps)
    Duration timestep = xiosHandler.getCalendarTimestep();
    metadata.setTime(xiosHandler.getCalendarStart());
    REQUIRE(xiosHandler.getCalendarStep() == 0);
    for (int ts = 1; ts <= 4; ts++) {

        // Update the current timestep and verify it's updated in XIOS
        metadata.incrementTime(timestep);
        REQUIRE(xiosHandler.getCalendarStep() == ts);

        // Update diagnostics
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {
                u(i, j) = ts;
            }
        }

        // Set up two ModelStates: one for restarts and one for diagnostics
        ModelState restarts = { {
                                    { maskName, mask },
                                    { coordsName, coordinates },
                                    { hiceName, hice },
                                    { ticeName, tice },
                                    { ciceName, cice },
                                },
            {} };
        ModelState diagnostics = { {
                                       { uName, u },
                                   },
            {} };

        // Write out diagnostics and then restarts
        pio->writeDiagnosticTime(diagnostics, diagnosticFilename);
        grid.dumpModelState(restarts, restartFilename, true);
    }

    // Check the files have indeed been created then remove it
    REQUIRE(std::filesystem::exists("xios_test_output_20230317171100-20230317201059.nc"));
    REQUIRE(std::filesystem::exists("xios_test_output_20230317201100-20230317231059.nc"));
    REQUIRE(std::filesystem::exists("xios_test_diagnostic.nc"));
    if (rank == 0) {
        std::filesystem::remove("xios_test_output_20230317171100-20230317201059.nc");
        std::filesystem::remove("xios_test_output_20230317201100-20230317231059.nc");
        std::filesystem::remove("xios_test_diagnostic.nc");
    }

    xiosHandler.context_finalize();
    Finalizer::finalize();
}

}
