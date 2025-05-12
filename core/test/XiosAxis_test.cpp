/*!
 * @file    XiosAxis_test.cpp
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @date    29 Apr 2025
 * @brief   Tests for XIOS axes
 * @details
 * This test is designed to test axis functionality of the C++ interface
 * for XIOS.
 *
 */
#include <doctest/extensions/doctest_mpi.h>
#undef INFO

#include "StructureModule/include/ParametricGrid.hpp"
#include "include/Finalizer.hpp"
#include "include/Xios.hpp"

namespace Nextsim {

/*!
 * TestXiosAxis
 *
 * This function tests the axis functionality of the C++ interface for XIOS. It
 * needs to be run with 2 ranks i.e.,
 *
 * `mpirun -n 2 ./testXiosAxis_MPI2`
 *
 */
MPI_TEST_CASE("TestXiosAxis", 2)
{
    enableXios();

    // Get the Xios singleton instance and check it's initialized
    Xios& xiosHandler = Xios::getInstance();
    REQUIRE(xiosHandler.isInitialized());
    REQUIRE(xiosHandler.getClientMPISize() == 2);

    // --- Tests for axis API
    const std::string axisId = { "axis_A" };
    REQUIRE_THROWS_WITH(xiosHandler.getAxisSize(axisId), "Xios: Undefined axis 'axis_A'");
    REQUIRE_THROWS_WITH(xiosHandler.getAxisValues(axisId), "Xios: Undefined axis 'axis_A'");
    xiosHandler.createAxis(axisId);
    REQUIRE_THROWS_WITH(xiosHandler.createAxis(axisId), "Xios: Axis 'axis_A' already exists");
    // Axis size
    REQUIRE_THROWS_WITH(xiosHandler.getAxisSize(axisId), "Xios: Undefined size for axis 'axis_A'");
    const size_t axisSize { 2 };
    xiosHandler.setAxisSize(axisId, axisSize);
    REQUIRE(xiosHandler.getAxisSize(axisId) == axisSize);
    // Axis values
    REQUIRE_THROWS_WITH(
        xiosHandler.getAxisValues(axisId), "Xios: Undefined values for axis 'axis_A'");
    REQUIRE_THROWS_WITH(xiosHandler.setAxisValues(axisId, { 0.0, 1.0, 2.0 }),
        "Xios: Size incompatible with values for axis 'axis_A'");
    std::vector<double> axisValues { 0.0, 1.0 };
    xiosHandler.setAxisValues(axisId, axisValues);
    std::vector<double> axis_A = xiosHandler.getAxisValues(axisId);
    REQUIRE(axis_A[0] == doctest::Approx(0.0));
    REQUIRE(axis_A[1] == doctest::Approx(1.0));

    xiosHandler.close_context_definition();
    xiosHandler.context_finalize();
    Finalizer::finalize();
}
}
