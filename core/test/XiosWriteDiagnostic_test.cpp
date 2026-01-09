/*!
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @brief   Tests for XIOS functionality for writing diagnostic files.
 * @details The functionality of writing diagnostics via XIOS is tested.
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
const std::string diagnosticFilename = testFilesDir + "/xios_test_diagnostic.nc";

static const int DGCOMP = 6;
static const int DGSTRESSCOMP = 8;
static const int CGDEGREE = 2;

namespace Nextsim {

/*!
 * TestXiosWriteDiagnostic
 *
 * Test writing of diagnostics via `writeDiagnosticTime`.
 */
MPI_TEST_CASE("TestXiosWriteDiagnostic", 2)
{
    std::stringstream config;
    config << "[model]" << std::endl;
    config << "start = 2023-03-17T17:11:00Z" << std::endl;
    config << "stop = 2023-03-17T23:11:00Z" << std::endl;
    config << "time_step = P0-0T01:30:00" << std::endl;
    config << "init_file = " << inputFilename << std::endl;
    config << "partition_file = xios_test_partition_metadata_2.nc" << std::endl;
    config << "restart_period = P0-0T03:00:00" << std::endl;
    config << "[XiosDiagnostic]" << std::endl;
    config << "filename = " << diagnosticFilename << std::endl;
    config << "field_names = " << hsnowName << std::endl;
    config << "period = P0-0T03:00:00" << std::endl;
    config << "split_period = P0-0T03:00:00" << std::endl;
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

    // Set ModelArray dimensions
    const size_t nx = ModelArray::size(ModelArray::Dimension::X);
    const size_t ny = ModelArray::size(ModelArray::Dimension::Y);
    REQUIRE(ModelArray::nComponents(ModelArray::Type::DG) == DGCOMP);
    REQUIRE(ModelArray::size(ModelArray::Dimension::DG) == DGCOMP);

    // Create ParametricGrid and ParaGridIO instances
    // NOTE: XIOS axes, domains, and grids are created by the ParaGridIO constructor
    Module::setImplementation<IStructure>("Nextsim::ParametricGrid");
    ParametricGrid grid;
    ParaGridIO* pio = new ParaGridIO(grid);
    grid.setIO(pio);

    // Set field type for diagnostics
    xiosHandler.setFieldType(hsnowName, ModelArray::Type::H);

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
    DGField hice(ModelArray::Type::DG);
    hice.resize();
    DGSField tice(ModelArray::Type::DGSTRESS);
    tice.resize();
    CGField uice(ModelArray::Type::CG);
    uice.resize();
    HField hsnow(ModelArray::Type::H);
    hsnow.resize();

    // Check files with the expected names don't exist yet
    REQUIRE_FALSE(std::filesystem::exists("xios_test_diagnostic*.nc"));

    // Check calendar step is zero initially
    REQUIRE(xiosHandler.getCalendarStep() == 0);

    // Simulate 4 iterations (timesteps)
    ModelMetadata& metadata = ModelMetadata::getInstance();
    const Duration& timestep = metadata.stepLength();
    int rank;
    MPI_Comm_rank(test_comm, &rank);
    for (int ts = 0; ts <= 4; ts++) {

        // Update diagnostics
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {
                hsnow(i, j) = 0.1 * ts;
            }
        }

        // Set up ModelStates for diagnostics and write out
        ModelState diagnostics = { {
                                       { hsnowName, hsnow },
                                   },
            {} };
        pio->writeDiagnosticTime(diagnostics, diagnosticFilename);

        // Update the current timestep and verify it's updated in XIOS
        metadata.incrementTime(timestep);
        REQUIRE(xiosHandler.getCalendarStep() == ts + 1);
    }

    // Check the files have indeed been created
    // NOTE: Don't remove them because their contents are checked in
    // XiosReadxios_test_diagnostic_test
    REQUIRE(std::filesystem::exists("xios_test_diagnostic_20230317171100-20230317201059.nc"));
    REQUIRE(std::filesystem::exists("xios_test_diagnostic_20230317201100-20230317231059.nc"));

    xiosHandler.context_finalize();
    Finalizer::finalize();
}

}
