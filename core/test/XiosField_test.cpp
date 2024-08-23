/*!
 * @file    XiosField_test.cpp
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @date    21 August 2024
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
    Duration timestep("P0-0T01:00:00");
    xios_handler.setCalendarTimestep(timestep);

    // Create an axis with two points
    xios_handler.createAxis("axis_A");
    xios_handler.setAxisValues("axis_A", { 0.0, 1.0 });

    // Create a 1D grid comprised of the single axis
    xios_handler.createGrid("grid_1D");
    xios_handler.gridAddAxis("grid_1D", "axis_A");

    // --- Tests for field API
    const std::string fieldId = "field_A";
    REQUIRE_THROWS_WITH(xios_handler.getFieldName(fieldId), "Xios: Undefined field 'field_A'");
    xios_handler.createField(fieldId);
    REQUIRE_THROWS_WITH(xios_handler.createField(fieldId), "Xios: Field 'field_A' already exists");
    // Field name
    REQUIRE_THROWS_WITH(
        xios_handler.getFieldName(fieldId), "Xios: Undefined name for field 'field_A'");
    const std::string fieldName = "test_field";
    xios_handler.setFieldName(fieldId, fieldName);
    REQUIRE(xios_handler.getFieldName(fieldId) == fieldName);
    // Operation
    REQUIRE_THROWS_WITH(
        xios_handler.getFieldOperation(fieldId), "Xios: Undefined operation for field 'field_A'");
    const std::string operation = "instant";
    xios_handler.setFieldOperation(fieldId, operation);
    REQUIRE(xios_handler.getFieldOperation(fieldId) == operation);
    // Grid reference
    REQUIRE_THROWS_WITH(xios_handler.getFieldGridRef(fieldId),
        "Xios: Undefined grid reference for field 'field_A'");
    const std::string gridRef = "grid_1D";
    xios_handler.setFieldGridRef(fieldId, gridRef);
    REQUIRE(xios_handler.getFieldGridRef(fieldId) == gridRef);
    // Read access
    const bool readAccess(true);
    xios_handler.setFieldReadAccess(fieldId, readAccess);
    REQUIRE(xios_handler.getFieldReadAccess(fieldId));
    // Frequency offset
    Duration freqOffset = timestep;
    xios_handler.setFieldFreqOffset(fieldId, freqOffset);
    // TODO: Overload == for Duration
    REQUIRE(xios_handler.getFieldFreqOffset(fieldId).seconds() == freqOffset.seconds());

    xios_handler.close_context_definition();
    xios_handler.context_finalize();
}
}
