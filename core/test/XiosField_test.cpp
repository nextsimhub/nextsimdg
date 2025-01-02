/*!
 * @file    XiosField_test.cpp
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @date    02 Jan 2025
 * @brief   Tests for XIOS axes
 * @details
 * This test is designed to test axis functionality of the C++ interface
 * for XIOS.
 *
 */
#include <doctest/extensions/doctest_mpi.h>
#undef INFO

#include "StructureModule/include/ParametricGrid.hpp"
#include "include/Xios.hpp"

namespace Nextsim {

/*!
 * TestXiosField
 *
 * This function tests the field functionality of the C++ interface for XIOS. It
 * needs to be run with 2 ranks i.e.,
 *
 * `mpirun -n 2 ./testXiosField_MPI2`
 *
 */
MPI_TEST_CASE("TestXiosField", 2)
{
    enableXios();

    // Get the Xios singleton instance and check it's initialized
    Xios* xiosHandler = Xios::getInstance();
    REQUIRE(xiosHandler->isInitialized());
    const size_t size = xiosHandler->getClientMPISize();
    REQUIRE(size == 2);
    const size_t rank = xiosHandler->getClientMPIRank();

    // Create an axis with two points
    xiosHandler->createAxis("axis_A");
    xiosHandler->setAxisValues("axis_A", { 0.0, 1.0 });

    // Create a 1D grid comprised of the single axis
    xiosHandler->createGrid("grid_1D");
    xiosHandler->gridAddAxis("grid_1D", "axis_A");

    // --- Tests for field API
    const std::string fieldId = "field_A";
    REQUIRE_THROWS_WITH(xiosHandler->getFieldName(fieldId), "Xios: Undefined field 'field_A'");
    xiosHandler->createField(fieldId);
    REQUIRE_THROWS_WITH(xiosHandler->createField(fieldId), "Xios: Field 'field_A' already exists");
    // Field name
    REQUIRE_THROWS_WITH(
        xiosHandler->getFieldName(fieldId), "Xios: Undefined name for field 'field_A'");
    const std::string fieldName = "test_field";
    xiosHandler->setFieldName(fieldId, fieldName);
    REQUIRE(xiosHandler->getFieldName(fieldId) == fieldName);
    // Operation
    REQUIRE_THROWS_WITH(
        xiosHandler->getFieldOperation(fieldId), "Xios: Undefined operation for field 'field_A'");
    const std::string operation = "instant";
    xiosHandler->setFieldOperation(fieldId, operation);
    REQUIRE(xiosHandler->getFieldOperation(fieldId) == operation);
    // Grid reference
    REQUIRE_THROWS_WITH(xiosHandler->getFieldGridRef(fieldId),
        "Xios: Undefined grid reference for field 'field_A'");
    const std::string gridRef = "grid_1D";
    xiosHandler->setFieldGridRef(fieldId, gridRef);
    REQUIRE(xiosHandler->getFieldGridRef(fieldId) == gridRef);
    // Read access
    const bool readAccess(true);
    xiosHandler->setFieldReadAccess(fieldId, readAccess);
    REQUIRE(xiosHandler->getFieldReadAccess(fieldId));
    // Frequency offset
    Duration freqOffset = xiosHandler->getCalendarTimestep();
    xiosHandler->setFieldFreqOffset(fieldId, freqOffset);
    // TODO: Overload == for Duration
    REQUIRE(xiosHandler->getFieldFreqOffset(fieldId).seconds() == freqOffset.seconds());

    xiosHandler->close_context_definition();
    xiosHandler->context_finalize();
}
}
