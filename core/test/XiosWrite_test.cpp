/*!
 * @file    XiosWrite_test.cpp
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @date    27 Aug 2024
 * @brief   Tests for XIOS write method
 * @details
 * This test is designed to test the write method of the C++ interface
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
 * TestXiosWrite
 *
 * This function tests the file writing functionality of the C++ interface for XIOS. It
 * needs to be run with 2 ranks i.e.,
 *
 * `mpirun -n 2 ./testXiosWrite_MPI2`
 *
 */
MPI_TEST_CASE("TestXiosWrite", 2)
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

    // Create fields on the two grids
    xios_handler.createField("field_2D");
    xios_handler.setFieldOperation("field_2D", "instant");
    xios_handler.setFieldGridRef("field_2D", "grid_2D");
    xios_handler.setFieldReadAccess("field_2D", false);
    xios_handler.createField("field_3D");
    xios_handler.setFieldOperation("field_3D", "instant");
    xios_handler.setFieldGridRef("field_3D", "grid_3D");
    xios_handler.setFieldReadAccess("field_3D", false);

    // Create an output file to hold data from both fields
    xios_handler.createFile("xios_test_output");
    xios_handler.setFileType("xios_test_output", "one_file");
    xios_handler.setFileOutputFreq("xios_test_output", timestep);
    xios_handler.setFileSplitFreq("xios_test_output", Duration("P0-0T03:00:00"));
    xios_handler.setFileMode("xios_test_output", "write");
    xios_handler.fileAddField("xios_test_output", "field_2D");
    xios_handler.fileAddField("xios_test_output", "field_3D");

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

    // --- Tests for writing to file
    // Create some fake data to test writing methods
    HField field_2D(ModelArray::Type::H);
    field_2D.resize();
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            field_2D(i, j) = 1.0 * (i + nx * j);
        }
    }
    HField field_3D(ModelArray::Type::Z);
    field_3D.resize();
    for (size_t k = 0; k < nz; ++k) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {
                field_3D(i, j, k) = 1.0 * (i + nx * (j + ny * k));
            }
        }
    }
    // Verify calendar step is starting from zero
    REQUIRE(xios_handler.getCalendarStep() == 0);
    // Check a file with the expected name doesn't exist yet
    REQUIRE_FALSE(std::filesystem::exists("xios_test_output*.nc"));
    // Simulate 4 iterations (timesteps)
    for (int ts = 1; ts <= 4; ts++) {
        // Update the current timestep
        xios_handler.updateCalendar(ts);
        // Send data to XIOS to be written to disk
        pio->write("field_2D", field_2D);
        pio->write("field_3D", field_3D);
        // Verify timestep
        REQUIRE(xios_handler.getCalendarStep() == ts);
    }
    // Check the files have indeed been created then remove it
    REQUIRE(std::filesystem::exists("xios_test_output_20230317171100-20230317201059.nc"));
    REQUIRE(std::filesystem::exists("xios_test_output_20230317201100-20230317231059.nc"));
    if (rank == 0) {
        std::filesystem::remove("xios_test_output_20230317171100-20230317201059.nc");
        std::filesystem::remove("xios_test_output_20230317201100-20230317231059.nc");
    }
    xios_handler.context_finalize();

    // TODO: Consider adding a 4D test case
}
}
