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

static const int DG = 3;

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
    config << "init_file = xios_test_input.nc" << std::endl;
    config << "restart_period = P0-0T01:30:00" << std::endl;
    config << "partition_file = xios_test_partition_metadata_2.nc" << std::endl;
    config << "[XiosInput]" << std::endl;
    config << "field_names = " << maskName << "," << coordsName << "," << hiceName << std::endl;
    config << "[XiosForcing]" << std::endl;
    config << "filename = xios_test_forcing.nc" << std::endl;
    config << "field_names = " << uName << std::endl;
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

    // Create ParametricGrid and ParaGridIO instances
    Module::setImplementation<IStructure>("Nextsim::ParametricGrid");
    ParametricGrid grid;
    ParaGridIO* pio = new ParaGridIO(grid);
    grid.setIO(pio);

    // Get the Xios singleton instance and check it's initialized
    Xios& xiosHandler = Xios::getInstance();
    REQUIRE(xiosHandler.isInitialized());
    REQUIRE(xiosHandler.getClientMPISize() == 2);

    // TODO: We could deduce this from the NetCDF file
    ModelArray::setNComponents(ModelArray::Type::DG, DG);
    ModelArray::setNComponents(ModelArray::Type::VERTEX, ModelArray::nCoords);
    REQUIRE(ModelArray::nComponents(ModelArray::Type::DG) == DG);
    REQUIRE(ModelArray::nComponents(ModelArray::Type::VERTEX) == ModelArray::nCoords);

    // Affix ModelMetadata to Xios handler
    // TODO: Automate this - can't be inlined in Xios::getInstance because need set field types
    xiosHandler.affixModelMetadata();

    xiosHandler.close_context_definition();

    // Check calendar step is zero initially
    REQUIRE(xiosHandler.getCalendarStep() == 0);

    // Check the input file exists
    REQUIRE(std::filesystem::exists(filename));

    // Deduce the local lengths of the two dimensions
    const size_t nx = ModelArray::size(ModelArray::Dimension::X);
    const size_t ny = ModelArray::size(ModelArray::Dimension::Y);

    // Read restarts from file and check they take the expected values
    ModelMetadata& metadata = ModelMetadata::getInstance();
    metadata.setTime(xiosHandler.getCalendarStart());
    REQUIRE(xiosHandler.getCalendarStep() == 0);
    ModelState restarts = grid.getModelState(filename);
    for (auto& entry : restarts.data) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {
                if (entry.first == maskName) {
                    REQUIRE(entry.second(i, j) == doctest::Approx(j >= 1 ? 1.0 : 0.0));
                } else if (entry.first == coordsName) {
                    REQUIRE(entry.second.components({ i, j })[0] == doctest::Approx(i));
                    REQUIRE(entry.second.components({ i, j })[1] == doctest::Approx(j));
                } else if (entry.first == hiceName) {
                    for (size_t d = 0; d < DG; ++d) {
                        float expected = 1.0 * (d + DG * (i + nx * j));
                        REQUIRE(entry.second.components({ i, j })[d] == doctest::Approx(expected));
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
        ModelState forcings = pio->readForcingTimeStatic(forcingFieldNames, time, filename);
        for (auto& entry : forcings.data) {
            for (size_t j = 0; j < ny; ++j) {
                for (size_t i = 0; i < nx; ++i) {
                    REQUIRE(entry.second(i, j) == doctest::Approx(ts));
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
