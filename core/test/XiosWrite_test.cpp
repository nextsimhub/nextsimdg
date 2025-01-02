/*!
 * @file    XiosWrite_test.cpp
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @date    02 Jan 2025
 * @brief   Tests for XIOS write functionality
 * @details
 * This test is designed to test the file writing functionality of the C++
 * interface for XIOS.
 *
 */
#include <doctest/extensions/doctest_mpi.h>
#undef INFO

#include "StructureModule/include/ParametricGrid.hpp"
#include "include/ModelMetadata.hpp"
#include "include/NextsimModule.hpp"
#include "include/ParaGridIO.hpp"
#include "include/Xios.hpp"
#include "include/gridNames.hpp"

#include <filesystem>

namespace Nextsim {

/*!
 * TestXiosWrite
 *
 * This function tests the file writing functionality of the C++ interface for XIOS for fields with
 * two and three spatial dimensions. The test runs with two MPI ranks.
 */
MPI_TEST_CASE("TestXiosWrite", 2)
{
    enableXios();

    // Create ParametricGrid and ParaGridIO instances
    Module::setImplementation<IStructure>("Nextsim::ParametricGrid");
    ParametricGrid grid;
    ParaGridIO* pio = new ParaGridIO(grid);
    grid.setIO(pio);

    // Create Xios singleton instance and check it's initialized
    REQUIRE(Xios::isNull());
    Xios* xiosHandler = Xios::getInstance("P0-0T01:30:00", "test", "2023-03-17T17:11:00Z");
    REQUIRE(xiosHandler->isInitialized());
    const size_t size = xiosHandler->getClientMPISize();
    REQUIRE(size == 2);
    const size_t rank = xiosHandler->getClientMPIRank();

    // Set ModelArray dimensions
    const size_t nx_glo = 4;
    const size_t ny_glo = 2;
    const size_t nx = 2;
    const size_t ny = 2;
    const size_t nz = 2;
    ModelArray::setDimension(ModelArray::Dimension::X, nx_glo, nx, 0);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny_glo, ny, 0);
    ModelArray::setDimension(ModelArray::Dimension::Z, nz, nz, 0);

    // Create some fake data to test writing methods
    HField hice(ModelArray::Type::H);
    hice.resize();
    HField cice(ModelArray::Type::Z);
    cice.resize();
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            hice(i, j) = 1.0 * (i + nx * j);
        }
    }
    for (size_t k = 0; k < nz; ++k) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {
                cice(i, j, k) = 1.0 * (i + nx * (j + ny * k));
            }
        }
    }

    // Setup ModelState with fields above
    ModelState state = { {
                             { hiceName, hice },
                             { ciceName, cice },
                         },
        {} };
    pio->setupXios(state, "xios_test_output", false);

    // Check calendar step is zero initially
    REQUIRE(xiosHandler->getCalendarStep() == 0);

    // Check a file with the expected name doesn't exist yet
    REQUIRE_FALSE(std::filesystem::exists("xios_test_output*.nc"));

    // Create ModelMetadata instance
    ModelMetadata metadata;
    metadata.setTime(xiosHandler->getCalendarStart());

    // Simulate 4 iterations (timesteps)
    Duration timestep = xiosHandler->getCalendarTimestep();
    for (int ts = 1; ts <= 4; ts++) {
        // Update the current timestep and verify it's updated in XIOS
        metadata.incrementTime(timestep);
        REQUIRE(xiosHandler->getCalendarStep() == ts);
        // Send data to XIOS to be written to disk
        grid.dumpModelState(state, metadata, "xios_test_output.nc", true);
    }

    // Check the files have indeed been created then remove it
    REQUIRE(std::filesystem::exists("xios_test_output_20230317171100-20230317201059.nc"));
    REQUIRE(std::filesystem::exists("xios_test_output_20230317201100-20230317231059.nc"));
    if (xiosHandler->getClientMPIRank() == 0) {
        std::filesystem::remove("xios_test_output_20230317171100-20230317201059.nc");
        std::filesystem::remove("xios_test_output_20230317201100-20230317231059.nc");
    }

    xiosHandler->context_finalize();
}

}
