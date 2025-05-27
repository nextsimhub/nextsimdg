/*!
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @author  Adeleke Bankole <ab3191@cam.ac.uk>
 * @brief   Tests for XIOS grid
 * @details
 * This test is designed to test grid functionality of the C++ interface
 * for XIOS.
 *
 */
#include <doctest/extensions/doctest_mpi.h>
#undef INFO

#include "include/Finalizer.hpp"
#include "include/ModelMetadata.hpp"
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
    enableXios();

    // Get the Xios singleton instance and check it's initialized
    Xios& xiosHandler = Xios::getInstance();
    REQUIRE(xiosHandler.isInitialized());
    const size_t size = xiosHandler.getClientMPISize();
    REQUIRE(size == 2);
    const size_t rank = xiosHandler.getClientMPIRank();

    // Create ModelMetadata instance, which will create the domain
    const std::string gridId = "grid_2D";
    REQUIRE_THROWS_WITH(xiosHandler.getGridAxisIds(gridId), "Xios: Undefined grid 'grid_2D'");
    ModelMetadata metadata("xios_test_partition_metadata_2.nc", test_comm);
    REQUIRE_THROWS_WITH(xiosHandler.createGrid(gridId), "Xios: Grid 'grid_2D' already exists");

    // Add a vertical axis, too
    const std::string axisId = "z_axis";
    xiosHandler.createAxis(axisId);
    xiosHandler.setAxisValues(axisId, { 0.0, 1.0 });
    xiosHandler.gridAddAxis("grid_2D", axisId);
    std::vector<std::string> axisIds = xiosHandler.getGridAxisIds(gridId);
    REQUIRE(axisIds.size() == 1);
    REQUIRE(axisIds[0] == axisId);

    xiosHandler.close_context_definition();
    xiosHandler.context_finalize();
    Finalizer::finalize();
}
}
