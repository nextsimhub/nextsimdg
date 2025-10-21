/*!
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @author  Adeleke Bankole <ab3191@cam.ac.uk>
 * @brief   Tests for XIOS grid
 * @details
 * This test is designed to test grid functionality of the C++ interface
 * for XIOS.
 */
#include <doctest/extensions/doctest_mpi.h>
#undef INFO

#include "include/Finalizer.hpp"
#include "include/Model.hpp"
#include "include/ModelMPI.hpp"
#include "include/Xios.hpp"

namespace Nextsim {

/*!
 * TestXiosGrid
 *
 * This function tests the grid functionality of the C++ interface for XIOS. It
 * needs to be run with 2 ranks i.e.,
 *
 * `mpirun -n 2 ./testXiosGrid_MPI2`
 *
 */
MPI_TEST_CASE("TestXiosGrid", 2)
{
    // Enable XIOS in the 'config' and provide parameters to configure it
    enableXios();
    std::stringstream config;
    config << "[model]" << std::endl;
    config << "start = 2023-03-17T17:11:00Z" << std::endl;
    config << "stop = 2023-03-17T18:11:00Z" << std::endl;
    config << "time_step = P0-0T01:00:00" << std::endl;
    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    // Create ModelMetadata instance based off a partition metadata file
    auto& modelMPI = ModelMPI::getInstance(test_comm);
    auto& metadata = ModelMetadata::getInstance("xios_test_partition_metadata_2.nc");

    // Create a Model and configure it so that time options are parsed
    Model model;
    model.configureTime();

    // Get the Xios singleton instance and check it's initialized
    Xios& xiosHandler = Xios::getInstance();
    REQUIRE(xiosHandler.isInitialized());

    // Affix ModelMetadata to Xios handler
    // TODO: Automate this - can't be inlined in Xios::getInstance because need set field types
    xiosHandler.affixModelMetadata();

    const std::string gridId = "HGrid";
    REQUIRE_THROWS_WITH(xiosHandler.createGrid(gridId), "Xios: Grid 'HGrid' already exists");

    // Add a vertical axis, too
    const std::string axisId = "z_axis";
    xiosHandler.createAxis(axisId);
    xiosHandler.setAxisValues(axisId, { 0.0, 1.0 });
    xiosHandler.gridAddAxis("HGrid", axisId);
    std::vector<std::string> axisIds = xiosHandler.getGridAxisIds(gridId);
    REQUIRE(axisIds.size() == 1);
    REQUIRE(axisIds[0] == axisId);

    xiosHandler.close_context_definition();
    xiosHandler.context_finalize();
    Finalizer::finalize();
}
}
