/*!
 * @file    XiosField_test.cpp
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @author  Adeleke Bankole <ab3191@cam.ac.uk>
 * @date    07 May 2025
 * @brief   Tests for XIOS fields
 * @details
 * This test is designed to test field functionality of the C++ interface
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
 * TestXiosField
 *
 * This function tests the field functionality of the C++ interface for XIOS. It
 * needs to be run with 3 ranks i.e.,
 *
 * `mpirun -n 3 ./testXiosField_MPI3`
 *
 */
MPI_TEST_CASE("TestXiosField", 3)
{
    // Enable XIOS in the 'config' and provide parameters to configure it
    enableXios();
    std::stringstream config;
    config << "[XiosInput]" << std::endl;
    config << "field_names = field_C,field_D" << std::endl;
    config << "[XiosOutput]" << std::endl;
    config << "field_names = field_A,field_C" << std::endl;
    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    // Get the Xios singleton instance and check it's initialized
    Xios& xiosHandler = Xios::getInstance();
    REQUIRE(xiosHandler.isInitialized());
    const size_t size = xiosHandler.getClientMPISize();
    REQUIRE(size == 3);
    const size_t rank = xiosHandler.getClientMPIRank();

    // Create an axis with two points
    xiosHandler.createAxis("axis_A");
    xiosHandler.setAxisValues("axis_A", { 0.0, 1.0 });

    // Create a 1D grid comprised of the single axis
    xiosHandler.createGrid("grid_1D");
    xiosHandler.gridAddAxis("grid_1D", "axis_A");

    // --- Tests for field API
    const std::string fieldId = "field_A";
    REQUIRE_THROWS_WITH(xiosHandler.getFieldName(fieldId), "Xios: Undefined field 'field_A'");
    xiosHandler.createField(fieldId);
    REQUIRE_THROWS_WITH(xiosHandler.createField(fieldId), "Xios: Field 'field_A' already exists");
    // Disallow creation of fields that aren't in either config section
    REQUIRE_THROWS_WITH(xiosHandler.createField("field_B"),
        "Xios: Field 'field_B' cannot be found in the XiosInput or XiosOutput config sections");
    // Disallow fields that're in both config sections (for now)
    REQUIRE_THROWS_WITH(xiosHandler.createField("field_C"),
        "Xios: Field 'field_C' found in both the XiosInput and XiosOutput config sections");
    // Field name
    REQUIRE_THROWS_WITH(
        xiosHandler.getFieldName(fieldId), "Xios: Undefined name for field 'field_A'");
    const std::string fieldName = "test_field";
    xiosHandler.setFieldName(fieldId, fieldName);
    REQUIRE(xiosHandler.getFieldName(fieldId) == fieldName);
    // Operation
    REQUIRE_THROWS_WITH(
        xiosHandler.getFieldOperation(fieldId), "Xios: Undefined operation for field 'field_A'");
    const std::string operation = "instant";
    xiosHandler.setFieldOperation(fieldId, operation);
    REQUIRE(xiosHandler.getFieldOperation(fieldId) == operation);
    // Grid reference
    REQUIRE_THROWS_WITH(
        xiosHandler.getFieldGridRef(fieldId), "Xios: Undefined grid reference for field 'field_A'");
    const std::string gridRef = "grid_1D";
    xiosHandler.setFieldGridRef(fieldId, gridRef);
    REQUIRE(xiosHandler.getFieldGridRef(fieldId) == gridRef);
    // Read access
    REQUIRE(!xiosHandler.getFieldReadAccess(fieldId));
    xiosHandler.createField("field_D");
    xiosHandler.setFieldGridRef("field_D", gridRef);
    REQUIRE(xiosHandler.getFieldReadAccess("field_D"));
    // Frequency offset
    Duration freqOffset = xiosHandler.getCalendarTimestep();
    xiosHandler.setFieldFreqOffset(fieldId, freqOffset);
    // TODO: Overload == for Duration
    REQUIRE(xiosHandler.getFieldFreqOffset(fieldId).seconds() == freqOffset.seconds());

    xiosHandler.close_context_definition();
    xiosHandler.context_finalize();
    Finalizer::finalize();
}
}
