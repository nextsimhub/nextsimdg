/*!
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @author  Adeleke Bankole <ab3191@cam.ac.uk>
 * @brief   Tests for XIOS axes
 * @details
 * This test is designed to test axis functionality of the C++ interface
 * for XIOS.
 */
#include <doctest/extensions/doctest_mpi.h>
#undef INFO

#include "StructureModule/include/ParametricGrid.hpp"
#include "include/Finalizer.hpp"
#include "include/Model.hpp"
#include "include/ModelMPI.hpp"
#include "include/Xios.hpp"

namespace Nextsim {

/*!
 * TestXiosAxis
 *
 * This function tests the axis functionality of the C++ interface for XIOS. It
 * needs to be run with 3 ranks i.e.,
 *
 * `mpirun -n 3 ./testXiosAxis_MPI3`
 *
 */
MPI_TEST_CASE("TestXiosAxis", 3)
{
    std::stringstream config;
    config << "[model]" << std::endl;
    config << "start = 2023-03-17T17:11:00Z" << std::endl;
    config << "stop = 2023-03-17T18:11:00Z" << std::endl;
    config << "time_step = P0-0T01:00:00" << std::endl;
    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    // Create ModelMetadata instance based off a partition metadata file
    auto& modelMPI = ModelMPI::getInstance(test_comm);
    auto& metadata = ModelMetadata::getInstance("xios_test_partition_metadata_3.nc");

    // Create a Model and configure it so that time options are parsed
    Model model;
    model.configureTime();

    // Get the Xios singleton instance and check it's initialized
    // NOTE: The singleton is created during configureTime
    Xios& xiosHandler = Xios::getInstance();
    REQUIRE(xiosHandler.isInitialized());

    // --- Tests for axis API
    const std::string axisId = { "axis_A" };
    REQUIRE_THROWS_WITH(xiosHandler.getAxisSize(axisId), "Xios: Undefined axis 'axis_A'");
    xiosHandler.createAxis(axisId);
    REQUIRE_THROWS_WITH(xiosHandler.createAxis(axisId), "Xios: Axis 'axis_A' already exists");
    // Axis size
    REQUIRE_THROWS_WITH(xiosHandler.getAxisSize(axisId), "Xios: Undefined size for axis 'axis_A'");
    const size_t axisSize { 2 };
    xiosHandler.setAxisSize(axisId, axisSize);
    REQUIRE(xiosHandler.getAxisSize(axisId) == axisSize);

    xiosHandler.close_context_definition();
    xiosHandler.context_finalize();
    Finalizer::finalize();
}
}
