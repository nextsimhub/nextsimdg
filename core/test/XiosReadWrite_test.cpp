/*!
 * @file    XiosReadWrite_test.cpp
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @date    26 Nov 2024
 * @brief   Tests for XIOS write method
 * @details
 * This test is designed to test the read and write methods of the C++
 * interface for XIOS.
 *
 */
#include <doctest/extensions/doctest_mpi.h>
#undef INFO

#include "StructureModule/include/ParametricGrid.hpp"
#include "include/NextsimModule.hpp"
#include "include/ParaGridIO_Xios.hpp"
#include "include/Xios.hpp"

#include <filesystem>

namespace Nextsim {

/*!
 * Create a formatted attribute ID.
 *
 * @param label Label to use for the attribute format
 * @param dim Spatial dimensions
 * @return "<label>_<dim>D"
 */
std::string formatId(const std::string label, const int dim)
{
    return (boost::format("%1$s_%2$1dD") % label % dim).str();
}

/*!
 * Set up the XIOS handler class for testing file reading and writing.
 *
 * The function assumes two MPI ranks.
 *
 * @param xios_handler Pointer to an Xios handler class
 * @param dim The number of spatial dimensions
 * @param read If true, set up for file reading test, otherwise for file writing test
 * @return Appropriately configured Xios handler class instance
 */
ModelState setupXiosHandler(Xios* xios_handler, int dim, bool read)
{
    if ((dim != 2) && (dim != 3)) {
        throw std::invalid_argument("Test only implemented for 2D and 3D cases");
    }

    // Create ParametricGrid and ParaGridIO instances
    Module::setImplementation<IStructure>(
        "Nextsim::ParametricGrid"); // TODO-JGW: Why can't I pull PIO out?
    ParametricGrid grid;
    ParaGridIO* pio = new ParaGridIO(grid);
    grid.setIO(pio);

    // Check configuration of Xios handler
    pio->xiosHandler = xios_handler;
    REQUIRE(pio->xiosHandler->isInitialized());
    const size_t size = pio->xiosHandler->getClientMPISize();
    REQUIRE(size == 2);
    const size_t rank = pio->xiosHandler->getClientMPIRank();

    // Calendar setup
    pio->xiosHandler->setCalendarOrigin(TimePoint("2020-01-23T00:08:15Z"));
    pio->xiosHandler->setCalendarStart(TimePoint("2023-03-17T17:11:00Z"));
    Duration timestep("P0-0T01:30:00");
    pio->xiosHandler->setCalendarTimestep(timestep);

    // Set ModelArray dimensions corresponding to a 4x2 horizontal domain with a partition halving
    // the x-extent and a vertical axis with 2 points
    const size_t nx_glo = 4;
    const size_t ny_glo = 2;
    const size_t nx = 2;
    const size_t ny = 2;
    ModelArray::setDimension(ModelArray::Dimension::X, nx_glo, nx, 0);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny_glo, ny, 0);
    if (dim == 3) {
        const size_t nz = 2;
        ModelArray::setDimension(ModelArray::Dimension::Z, nz, nz, 0);
    }

    // TODO-JGW: setupXios should be called from dumpModelState or getModelState
    const std::string fieldId = formatId("field", dim);
    ModelState state;
    if (dim == 2) {
        HField field_2D(ModelArray::Type::H);
        field_2D.resize();
        state = { {
                      { fieldId, field_2D },
                  },
            {} };
    } else {
        HField field_3D(ModelArray::Type::Z);
        field_3D.resize();
        state = { {
                      { fieldId, field_3D },
                  },
            {} };
    }
    ModelMetadata meta; // TODO-JGW: Set up properly
    std::string fileId;
    if (read) {
        fileId = "xios_test_input";
    } else {
        fileId = "xios_test_output";
    }
    const std::string filename = fileId + ".nc";
    pio->setupXios(state, meta, filename, read);
    return state;
}

/*!
 * Test file reading for the Xios handler configuration in setupXiosHandler.
 *
 * @param xios_handler Pointer to Xios handler class instance configured using setupXiosHandler
 * @param field_A Reference to nextSIM-DG HField instance to test reading from file
 */
void readFile(Xios* xios_handler, HField& field_A, const std::string fieldId)
{
    // Verify calendar step is starting from zero
    REQUIRE(xios_handler->getCalendarStep() == 0);

    // Check the input file exists
    REQUIRE(std::filesystem::exists("xios_test_input.nc"));

    // Simulate 4 iterations (timesteps)
    for (int ts = 1; ts <= 4; ts++) {
        // Update the current timestep
        xios_handler->updateCalendar(ts);
        // Receive data from XIOS that is read from disk
        xios_handler->read(fieldId, field_A);
        // Verify timestep
        REQUIRE(xios_handler->getCalendarStep() == ts);
    }
}

/*!
 * Utility for checking that two double values are approximately equal.
 *
 * Without this (i.e., if it's inlined below) the first test passes but the second one fails. The
 * same is true with any REQUIRE call.
 *
 * @param val1 the first double
 * @param val2 the second double
 */
void assertIsClose(double val1, double val2) { REQUIRE(val1 == doctest::Approx(val2)); }

/*!
 * TestXiosRead_2D
 *
 * This function tests the file reading functionality of the C++ interface for XIOS for fields with
 * two spatial dimensions. The test runs with two MPI ranks.
 */
MPI_TEST_CASE("TestXiosRead_2D", 2)
{
    // TODO: Create XIOS handler along with ParaGridIO instance
    const std::string contextId = "read_2D";
    const std::string fieldId = "field_2D";
    Xios xios_handler(contextId);

    ModelState state = setupXiosHandler(&xios_handler, 2, true);
    ModelArray field_2D = state.data[fieldId];

    // Verify fields are read in correctly
    readFile(&xios_handler, field_2D, fieldId);
    const size_t nx = xios_handler.getDomainLocalXSize("xy_domain");
    const size_t ny = xios_handler.getDomainLocalYSize("xy_domain");
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            assertIsClose(field_2D(i, j), i + nx * j);
        }
    }
    xios_handler.context_finalize();
}

/*!
 * TestXiosRead_3D
 *
 * This function tests the file reading functionality of the C++ interface for XIOS for fields with
 * three spatial dimensions. The test runs with two MPI ranks.
 */
MPI_TEST_CASE("TestXiosRead_3D", 2)
{
    // TODO: Create XIOS handler along with ParaGridIO instance
    const std::string contextId = "read_3D";
    const std::string fieldId = "field_3D";
    Xios xios_handler(contextId);

    ModelState state = setupXiosHandler(&xios_handler, 3, true);
    ModelArray field_3D = state.data[fieldId];

    // Verify fields are read in correctly
    readFile(&xios_handler, field_3D, fieldId);
    const size_t nx = xios_handler.getDomainLocalXSize("xy_domain");
    const size_t ny = xios_handler.getDomainLocalYSize("xy_domain");
    const size_t nz = xios_handler.getAxisSize("z_axis");
    for (size_t k = 0; k < nz; ++k) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {
                assertIsClose(field_3D(i, j, k), i + nx * (j + ny * k));
            }
        }
    }
    xios_handler.context_finalize();
}

/*!
 * Test file writing for the Xios handler configuration in setupXiosHandler.
 *
 * @param xios_handler Pointer to Xios handler class instance configured using setupXiosHandler
 * @param state Reference to ModelState instance containing fields to be written to file
 */
void testFileWrite(Xios* xios_handler, ModelState& state)
{
    // Verify calendar step is starting from zero
    REQUIRE(xios_handler->getCalendarStep() == 0);

    // Check a file with the expected name doesn't exist yet
    REQUIRE_FALSE(std::filesystem::exists("xios_test_output*.nc"));

    // Simulate 4 iterations (timesteps)
    for (int ts = 1; ts <= 4; ts++) {
        // Update the current timestep
        xios_handler->updateCalendar(ts);
        // Send data to XIOS to be written to disk
        for (auto& entry : state.data) {
            xios_handler->write(entry.first, entry.second);
        }
        // Verify timestep
        REQUIRE(xios_handler->getCalendarStep() == ts);
    }

    // Check the files have indeed been created then remove it
    REQUIRE(std::filesystem::exists("xios_test_output_20230317171100-20230317201059.nc"));
    REQUIRE(std::filesystem::exists("xios_test_output_20230317201100-20230317231059.nc"));
    if (xios_handler->getClientMPIRank() == 0) {
        std::filesystem::remove("xios_test_output_20230317171100-20230317201059.nc");
        std::filesystem::remove("xios_test_output_20230317201100-20230317231059.nc");
    }

    xios_handler->context_finalize();
}

/*!
 * TestXiosWrite_2D
 *
 * This function tests the file writing functionality of the C++ interface for XIOS for fields with
 * two spatial dimensions. The test runs with two MPI ranks.
 */
MPI_TEST_CASE("TestXiosWrite_2D", 2)
{
    // TODO: Create XIOS handler along with ParaGridIO instance
    const std::string contextId = "write_2D";
    const std::string fieldId = "field_2D";
    Xios xios_handler(contextId);

    ModelState state = setupXiosHandler(&xios_handler, 2, false);
    ModelArray field_2D = state.data[fieldId];

    // Create some fake data to test writing methods
    const size_t nx = xios_handler.getDomainLocalXSize("xy_domain");
    const size_t ny = xios_handler.getDomainLocalYSize("xy_domain");
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            field_2D(i, j) = 1.0 * (i + nx * j);
        }
    }

    testFileWrite(&xios_handler, state);
}

/*!
 * TestXiosWrite_3D
 *
 * This function tests the file writing functionality of the C++ interface for XIOS for fields with
 * three spatial dimensions. The test runs with two MPI ranks.
 */
MPI_TEST_CASE("TestXiosWrite_3D", 2)
{
    // TODO: Create XIOS handler along with ParaGridIO instance
    const std::string contextId = "write_3D";
    const std::string fieldId = "field_3D";
    Xios xios_handler(contextId);

    ModelState state = setupXiosHandler(&xios_handler, 3, false);
    ModelArray field_3D = state.data[fieldId];

    // Create some fake data to test writing methods
    const size_t nx = xios_handler.getDomainLocalXSize("xy_domain");
    const size_t ny = xios_handler.getDomainLocalYSize("xy_domain");
    const size_t nz = xios_handler.getAxisSize("z_axis");
    for (size_t k = 0; k < nz; ++k) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {
                field_3D(i, j, k) = 1.0 * (i + nx * (j + ny * k));
            }
        }
    }

    testFileWrite(&xios_handler, state);
}

// TODO: Consider adding 4D test cases

}
