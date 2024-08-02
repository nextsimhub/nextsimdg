/*!
 * @file    XiosWrite_test.cpp
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @date    21 August 2024
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
    xios_handler.setCalendarTimestep(Duration("P0-0T01:30:00"));

    // Domain setup
    xios_handler.createDomain("xy_domain");
    xios_handler.setDomainType("xy_domain", "rectilinear");
    const size_t nx_glo { 4 };
    xios_handler.setDomainGlobalXSize("xy_domain", nx_glo);
    const size_t ny_glo { 2 };
    xios_handler.setDomainGlobalYSize("xy_domain", ny_glo);
    const size_t nx { nx_glo / size };
    xios_handler.setDomainLocalXSize("xy_domain", nx);
    const size_t ny { ny_glo };
    xios_handler.setDomainLocalYSize("xy_domain", ny);
    xios_handler.setDomainLocalXStart("xy_domain", nx * rank);
    xios_handler.setDomainLocalYStart("xy_domain", 0);
    std::vector<double> vx { -1.0 + rank, -0.5 + rank };
    xios_handler.setDomainLocalXValues("xy_domain", vx);
    std::vector<double> vy { -1.0, 1.0 };
    xios_handler.setDomainLocalYValues("xy_domain", vy);

    // Axis setup
    const int nz = 2;
    xios_handler.createAxis("z_axis");
    xios_handler.setAxisSize("z_axis", nz);
    xios_handler.setAxisValues("z_axis", { 0.0, 1.0 });

    // Grid setup
    xios_handler.createGrid("grid_2D");
    xios_handler.gridAddDomain("grid_2D", "xy_domain");
    xios_handler.createGrid("grid_3D");
    xios_handler.gridAddDomain("grid_3D", "xy_domain");
    xios_handler.gridAddAxis("grid_3D", "z_axis");

    // Field setup
    xios_handler.createField("field_2D");
    xios_handler.setFieldOperation("field_2D", "instant");
    xios_handler.setFieldGridRef("field_2D", "grid_2D");
    xios_handler.setFieldReadAccess("field_2D", false);
    xios_handler.createField("field_3D");
    xios_handler.setFieldOperation("field_3D", "instant");
    xios_handler.setFieldGridRef("field_3D", "grid_3D");
    xios_handler.setFieldReadAccess("field_3D", false);

    // File setup
    xios_handler.createFile("xios_test_output");
    xios_handler.setFileType("xios_test_output", "one_file");
    xios_handler.setFileOutputFreq("xios_test_output", Duration("P0-0T01:30:00"));
    xios_handler.setFileSplitFreq("xios_test_output", Duration("P0-0T03:00:00"));
    xios_handler.setFileMode("xios_test_output", "write");
    xios_handler.fileAddField("xios_test_output", "field_2D");
    xios_handler.fileAddField("xios_test_output", "field_3D");

    xios_handler.close_context_definition();

    // --- Tests for writing to file
    Module::setImplementation<IStructure>("Nextsim::ParametricGrid");
    ModelArray::setDimension(ModelArray::Dimension::X, nx, nx, 0);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny, ny, 0);
    ModelArray::setDimension(ModelArray::Dimension::Z, nz, nz, 0);
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
        xios_handler.write("field_2D", field_2D);
        xios_handler.write("field_3D", field_3D);
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
