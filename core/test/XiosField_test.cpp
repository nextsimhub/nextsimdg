/*!
 * @file    XiosField_test.cpp
 * @author  Joe Wallwork <jw2423@cam.ac.uk
 * @date    27 June 2024
 * @brief   Tests for XIOS axes
 * @details
 * This test is designed to test axis functionality of the C++ interface
 * for XIOS.
 *
 */
#include <doctest/extensions/doctest_mpi.h>
#undef INFO

#include "StructureModule/include/ParametricGrid.hpp"
#include "include/Configurator.hpp"
#include "include/Xios.hpp"

#include <iostream>

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

    // Enable XIOS in the 'config'
    Configurator::clearStreams();
    std::stringstream config;
    config << "[xios]" << std::endl << "enable = true" << std::endl;
    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    // Initialize an Xios instance called xios_handler
    Xios xios_handler;
    REQUIRE(xios_handler.isInitialized());
    const size_t size = xios_handler.getClientMPISize();
    REQUIRE(size == 2);
    const size_t rank = xios_handler.getClientMPIRank();

    // Set timestep as a minimum
    xios_handler.setCalendarTimestep(Duration("P0-0T01:00:00"));

    // Axis setup
    xios_handler.createAxis("axis_A");
    xios_handler.setAxisValues("axis_A", { 0, 1 });

    // Grid setup
    xios_handler.createGrid("grid_1D");
    xios_handler.gridAddAxis("grid_1D", "axis_A");

    // --- Tests for field API
    const std::string fieldId = { "field_A" };
    xios_handler.createField(fieldId);
    // Field name
    const std::string fieldName = { "test_field" };
    REQUIRE_FALSE(xios_handler.isDefinedFieldName(fieldId));
    xios_handler.setFieldName(fieldId, fieldName);
    REQUIRE(xios_handler.getFieldName(fieldId) == fieldName);
    REQUIRE(xios_handler.isDefinedFieldName(fieldId));
    // Operation
    const std::string operation = { "instant" };
    REQUIRE_FALSE(xios_handler.isDefinedFieldOperation(fieldId));
    xios_handler.setFieldOperation(fieldId, operation);
    REQUIRE(xios_handler.isDefinedFieldOperation(fieldId));
    REQUIRE(xios_handler.getFieldOperation(fieldId) == operation);
    // Grid reference
    const std::string gridRef = { "grid_2D" };
    REQUIRE_FALSE(xios_handler.isDefinedFieldGridRef(fieldId));
    xios_handler.setFieldGridRef(fieldId, gridRef);
    REQUIRE(xios_handler.isDefinedFieldGridRef(fieldId));
    REQUIRE(xios_handler.getFieldGridRef(fieldId) == gridRef);

    xios_handler.close_context_definition();
    xios_handler.context_finalize();
}
}
