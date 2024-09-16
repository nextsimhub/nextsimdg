/*!
 * @file    XiosRead_test.cpp
 * @author  Joe Wallwork <jw2423@cam.ac.uk
 * @date    27 Aug 2024
 * @brief   Tests for XIOS read method
 * @details
 * This test is designed to test the read method of the C++ interface
 * for XIOS.
 *
 */
#include <doctest/extensions/doctest_mpi.h>
#undef INFO

#include "StructureModule/include/ParametricGrid.hpp"
#include "include/Configurator.hpp"
#include "include/Module.hpp"
#include "include/ParaGridIO_Xios.hpp"
#include "include/Xios.hpp"

#include <filesystem>
#include <iostream>

namespace Nextsim {

/*!
 * TestXiosRead
 *
 * This function tests the file reading functionality of the C++ interface for XIOS. It
 * needs to be run with 2 ranks i.e.,
 *
 * `mpirun -n 2 ./testXiosRead_MPI2`
 *
 */
MPI_TEST_CASE("TestXiosRead", 2)
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

    // Calendar setup
    xios_handler.setCalendarOrigin(TimePoint("2020-01-23T00:08:15Z"));
    xios_handler.setCalendarStart(TimePoint("2023-03-17T17:11:00Z"));
    Duration timestep("P0-0T01:30:00");
    xios_handler.setCalendarTimestep(timestep);

    // Create a 4x2 horizontal domain with a partition halving the x-extent
    xios_handler.createDomain("xy_domain");
    xios_handler.setDomainType("xy_domain", "rectilinear");
    xios_handler.setDomainGlobalXSize("xy_domain", 4);
    xios_handler.setDomainGlobalYSize("xy_domain", 2);
    xios_handler.setDomainLocalXStart("xy_domain", 2 * rank);
    xios_handler.setDomainLocalYStart("xy_domain", 0);
    xios_handler.setDomainLocalXValues("xy_domain", { -1.0 + rank, -0.5 + rank });
    xios_handler.setDomainLocalYValues("xy_domain", { -1.0, 1.0 });

    // Create a vertical axis with 2 points
    xios_handler.createAxis("z_axis");
    xios_handler.setAxisValues("z_axis", { 0.0, 1.0 });

    // Create a 2D grid comprised of the xy-domain and a 3D grid which also includes the z-axis
    xios_handler.createGrid("grid_2D");
    xios_handler.gridAddDomain("grid_2D", "xy_domain");
    xios_handler.createGrid("grid_3D");
    xios_handler.gridAddDomain("grid_3D", "xy_domain");
    xios_handler.gridAddAxis("grid_3D", "z_axis");

    // Field setup
    xios_handler.createField("field_2D");
    xios_handler.setFieldOperation("field_2D", "instant");
    xios_handler.setFieldGridRef("field_2D", "grid_2D");
    xios_handler.setFieldReadAccess("field_2D", true);
    xios_handler.setFieldFreqOffset("field_2D", timestep);
    xios_handler.createField("field_3D");
    xios_handler.setFieldOperation("field_3D", "instant");
    xios_handler.setFieldGridRef("field_3D", "grid_3D");
    xios_handler.setFieldReadAccess("field_3D", true);
    xios_handler.setFieldFreqOffset("field_3D", timestep);

    // File setup
    xios_handler.createFile("xios_test_input");
    xios_handler.setFileType("xios_test_input", "one_file");
    xios_handler.setFileOutputFreq("xios_test_input", timestep);
    xios_handler.setFileMode("xios_test_input", "read");
    xios_handler.setFileParAccess("xios_test_input", "collective");
    xios_handler.fileAddField("xios_test_input", "field_2D");
    xios_handler.fileAddField("xios_test_input", "field_3D");

    // Set ModelArray dimensions
    Module::setImplementation<IStructure>("Nextsim::ParametricGrid");
    const size_t nx = xios_handler.getDomainLocalXSize("xy_domain");
    const size_t ny = xios_handler.getDomainLocalYSize("xy_domain");
    const size_t nz = xios_handler.getAxisSize("z_axis");
    ModelArray::setDimension(
        ModelArray::Dimension::X, xios_handler.getDomainGlobalXSize("xy_domain"), nx, 0);
    ModelArray::setDimension(
        ModelArray::Dimension::Y, xios_handler.getDomainGlobalYSize("xy_domain"), ny, 0);
    ModelArray::setDimension(ModelArray::Dimension::Z, nz, nz, 0);

    // Create ParametricGrid and ParaGridIO instances
    ParametricGrid grid;
    ParaGridIO* pio = new ParaGridIO(grid);
    grid.setIO(pio);

    xios_handler.close_context_definition();

    // --- Tests for reading to file
    // Create some fake data to test writing methods
    HField field_2D(ModelArray::Type::H);
    field_2D.resize();
    HField field_3D(ModelArray::Type::Z);
    field_3D.resize();
    // Verify calendar step is starting from zero
    REQUIRE(xios_handler.getCalendarStep() == 0);
    // Check the input file exists
    REQUIRE(std::filesystem::exists("xios_test_input.nc"));
    // Simulate 4 iterations (timesteps)
    for (int ts = 1; ts <= 4; ts++) {
        // Update the current timestep
        xios_handler.updateCalendar(ts);
        // Receive data from XIOS that is read from disk
        pio->read("field_2D", field_2D);
        pio->read("field_3D", field_3D);
        // Verify timestep
        REQUIRE(xios_handler.getCalendarStep() == ts);
    }
    // Verify fields have been read in correctly
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            REQUIRE(field_2D(i, j) == doctest::Approx(1.0 * (i + nx * j)));
        }
    }
    for (size_t k = 0; k < nz; ++k) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {
                REQUIRE(field_3D(i, j, k) == doctest::Approx(1.0 * (i + nx * (j + ny * k))));
            }
        }
    }
    xios_handler.context_finalize();

    // TODO: Consider adding a 4D test case
}
}
