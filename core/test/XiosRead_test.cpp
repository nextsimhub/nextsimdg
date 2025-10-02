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
const std::string filename = testFilesDir + "/xios_test_input.nc";

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
    config << "[XiosInput]" << std::endl;
    config << "filename = xios_test_input.nc" << std::endl;
    config << "field_names = hice" << std::endl;
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

    // TODO: Automate this
    xiosHandler.affixModelMetadata();

    // Create fields on the two grids
    // NOTE: Fields are created when the XIOS handler is constructed
    // NOTE: The 2D grid is created along with the 2D domain
    xiosHandler.setFieldOperation(hiceName, "instant");
    xiosHandler.setFieldGridRef(hiceName, "grid_2D");
    Duration timestep = xiosHandler.getCalendarTimestep();
    xiosHandler.setFieldFreqOffset(hiceName, timestep);

    xiosHandler.close_context_definition();

    // Create HField and ZField instances to read the data into
    HField hice(ModelArray::Type::H);
    hice.resize();

    // Setup ModelState with field above
    ModelState state = { {
                             { hiceName, hice },
                         },
        {} };

    // Check calendar step is zero initially
    REQUIRE(xiosHandler.getCalendarStep() == 0);

    // Check the input file exists
    REQUIRE(std::filesystem::exists(filename));

    // Deduce the local lengths of the two dimensions
    const size_t nx = ModelArray::definedDimensions.at(ModelArray::Dimension::X).localLength;
    const size_t ny = ModelArray::definedDimensions.at(ModelArray::Dimension::Y).localLength;

    // Simulate 4 iterations (timesteps)
    metadata.setTime(xiosHandler.getCalendarStart());
    for (int ts = 1; ts <= 4; ts++) {
        // Update the current timestep and verify it's updated in XIOS
        metadata.incrementTime(timestep);
        REQUIRE(xiosHandler.getCalendarStep() == ts);
        ModelState state = grid.getModelState(filename);
        for (auto& entry : state.data) {
            for (size_t j = 0; j < ny; ++j) {
                for (size_t i = 0; i < nx; ++i) {
                    REQUIRE(entry.second(i, j) == doctest::Approx(i + nx * j));
                }
            }
        }
    }
    xiosHandler.context_finalize();
    Finalizer::finalize();
}

}
