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
    config << "start = 2023-03-17T00:00:00Z" << std::endl;
    config << "stop = 2023-03-17T00:00:00Z" << std::endl;
    config << "time_step = P0-0T01:00:00" << std::endl;
    config << "partition_file = xios_test_partition_metadata_3.nc" << std::endl;
    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    // Create ModelMPI instance based off the test communicator
    auto& modelMPI = ModelMPI::getInstance(test_comm);

    // Create a Model and configure it so that time options are parsed
    Model model;
    model.configureTime();

    // Get the Xios singleton instance and check it's initialized
    // NOTE: The singleton is created when Xios::getInstance() is first called. In this test, this
    //       happens when the time sets set by ModelMetadata::setTime(). This occurs in the call to
    //       Model::configureTime() above.
    Xios& xiosHandler = Xios::getInstance();

    // Tests for axis size getter
    REQUIRE_THROWS_WITH(xiosHandler.getAxisSize("axis_A"), "Xios: Undefined axis 'axis_A'");
    REQUIRE(xiosHandler.getAxisSize("DGAxis") == 6);
    REQUIRE(xiosHandler.getAxisSize("DGSAxis") == 8);
    REQUIRE(xiosHandler.getAxisSize("VertexAxis") == 2);

    xiosHandler.context_finalize();
    Finalizer::finalize();
}
}
