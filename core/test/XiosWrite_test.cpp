/*!
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @brief   Tests for XIOS write functionality
 * @details
 * This test is designed to test the file writing functionality of the C++
 * interface for XIOS.
 */
#include <doctest/extensions/doctest_mpi.h>
#undef INFO

#include "StructureModule/include/ParametricGrid.hpp"
#include "include/Finalizer.hpp"
#include "include/ModelMetadata.hpp"
#include "include/NextsimModule.hpp"
#include "include/ParaGridIO.hpp"
#include "include/Xios.hpp"
#include "include/gridNames.hpp"

#include <filesystem>

static const int DG = 3;

namespace Nextsim {

/*!
 * TestXiosWrite
 *
 * This function tests the file writing functionality of the C++ interface for XIOS for fields with
 * two and three spatial dimensions. The test runs with two MPI ranks.
 */
MPI_TEST_CASE("TestXiosWrite", 2)
{
    // Enable XIOS in the 'config' and provide parameters to configure it
    enableXios();
    std::stringstream config;
    config << "[model]" << std::endl;
    config << "start = 2023-03-17T17:11:00Z" << std::endl;
    config << "time_step = P0-0T01:30:00" << std::endl;
    config << "[XiosOutput]" << std::endl;
    config << "period = P0-0T01:30:00" << std::endl;
    config << "filename = xios_test_output.nc" << std::endl;
    config << "field_names = " << maskName << "," << coordsName << "," << hiceName << std::endl;
    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    // Create ParametricGrid and ParaGridIO instances
    Module::setImplementation<IStructure>("Nextsim::ParametricGrid");
    ParametricGrid grid;
    ParaGridIO* pio = new ParaGridIO(grid);
    grid.setIO(pio);

    // Get the Xios singleton instance and check it's initialized
    Xios& xiosHandler = Xios::getInstance();
    REQUIRE(xiosHandler.isInitialized());
    REQUIRE(xiosHandler.getClientMPISize() == 2);

    // Set ModelArray dimensions
    const size_t nx_glo = 4;
    const size_t ny_glo = 2;
    const size_t nx = 2;
    const size_t ny = 2;
    ModelArray::setDimension(ModelArray::Dimension::X, nx_glo, nx, 0);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny_glo, ny, 0);
    ModelArray::setDimension(ModelArray::Dimension::XVERTEX, nx_glo + 1, nx + 1, 0);
    ModelArray::setDimension(ModelArray::Dimension::YVERTEX, ny_glo + 1, ny + 1, 0);
    ModelArray::setNComponents(ModelArray::Type::DG, DG);
    ModelArray::setNComponents(ModelArray::Type::VERTEX, ModelArray::nCoords);
    REQUIRE(ModelArray::nComponents(ModelArray::Type::DG) == DG);
    REQUIRE(ModelArray::nComponents(ModelArray::Type::VERTEX) == ModelArray::nCoords);

    // Create ModelMetadata instance based off a partition metadata file
    ModelMetadata metadata("xios_test_partition_metadata_2.nc", test_comm);
    xiosHandler.affixModelMetadata(metadata);

    // Create fields on the grid
    // NOTE: Fields are created when the XIOS handler is constructed
    // NOTE: The 2D grid is created along with the 2D domain
    Duration timestep = xiosHandler.getCalendarTimestep();
    for (std::string fieldName : { maskName, coordsName, hiceName }) {
        xiosHandler.setFieldOperation(fieldName, "instant");
        xiosHandler.setFieldFreqOffset(fieldName, timestep);
    }
    xiosHandler.setFieldType(maskName, ModelArray::Type::H);
    xiosHandler.setFieldType(coordsName, ModelArray::Type::VERTEX);
    xiosHandler.setFieldType(hiceName, ModelArray::Type::DG);

    // Set file split frequency
    // NOTE: Files are created when the XIOS handler is constructed
    const std::string fileId = "xios_test_output";
    xiosHandler.setFileSplitFreq(fileId, Duration("P0-0T03:00:00"));

    xiosHandler.close_context_definition();

    // Create some fake data to test writing methods
    HField mask(ModelArray::Type::H);
    mask.resize();
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            mask(i, j) = j >= 1 ? 1.0 : 0.0;
        }
    }
    VertexField coordinates(ModelArray::Type::VERTEX);
    coordinates.resize();
    for (size_t j = 0; j < ny + 1; ++j) {
        for (size_t i = 0; i < nx + 1; ++i) {
            coordinates.components({ i, j })[0] = (double)i;
            coordinates.components({ i, j })[1] = (double)j;
        }
    }
    DGField hice(ModelArray::Type::DG);
    hice.resize();
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            for (size_t d = 0; d < DG; ++d) {
                hice.components({ i, j })[d] = 1.0 * (d + DG * (i + nx * j));
            }
        }
    }

    // Setup ModelState with field above
    ModelState state = { {
                             { maskName, mask },
                             { coordsName, coordinates },
                             { hiceName, hice },
                         },
        {} };

    // Check calendar step is zero initially
    REQUIRE(xiosHandler.getCalendarStep() == 0);

    // Check a file with the expected name doesn't exist yet
    REQUIRE_FALSE(std::filesystem::exists("xios_test_output*.nc"));

    // Simulate 4 iterations (timesteps)
    metadata.setTime(xiosHandler.getCalendarStart());
    for (int ts = 1; ts <= 4; ts++) {
        // Update the current timestep and verify it's updated in XIOS
        metadata.incrementTime(timestep);
        REQUIRE(xiosHandler.getCalendarStep() == ts);
        grid.dumpModelState(state, metadata, "xios_test_output", true);
    }

    // Check the files have indeed been created then remove it
    REQUIRE(std::filesystem::exists("xios_test_output_20230317171100-20230317201059.nc"));
    REQUIRE(std::filesystem::exists("xios_test_output_20230317201100-20230317231059.nc"));
    if (xiosHandler.getClientMPIRank() == 0) {
        std::filesystem::remove("xios_test_output_20230317171100-20230317201059.nc");
        std::filesystem::remove("xios_test_output_20230317201100-20230317231059.nc");
    }

    xiosHandler.context_finalize();
    Finalizer::finalize();
}

}
