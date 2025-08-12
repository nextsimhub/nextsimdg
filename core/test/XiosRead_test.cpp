/*!
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @brief   Tests for XIOS read functionality
 * @details
 * This test is designed to test the file reading functionality of the C++
 * interface for XIOS.
 *
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

const std::string testFilesDir = TEST_FILES_DIR;
const std::string filename = testFilesDir + "/xios_test_input.nc";

static const int DG = 3;

namespace Nextsim {

/*!
 * TestXiosRead
 *
 * This function tests the file reading functionality of the C++ interface for XIOS for fields with
 * two and three spatial dimensions. The test runs with two MPI ranks.
 */
MPI_TEST_CASE("TestXiosRead", 2)
{
    // Enable XIOS in the 'config' and provide parameters to configure it
    enableXios();
    std::stringstream config;
    config << "[model]" << std::endl;
    config << "start = 2023-03-17T17:11:00Z" << std::endl;
    config << "time_step = P0-0T01:30:00" << std::endl;
    config << "[XiosInput]" << std::endl;
    config << "period = P0-0T01:30:00" << std::endl;
    config << "filename = xios_test_input.nc" << std::endl;
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

    // Create ModelMetadata instance based off a partition metadata file
    // NOTE: ModelArray dimensions are determined from the input file, if present
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
    // TODO: Automate the following
    xiosHandler.setFieldGridRef(maskName, "HGrid");
    xiosHandler.setFieldGridRef(hiceName, "DGGrid");
    xiosHandler.setFieldGridRef(coordsName, "VertexGrid");

    xiosHandler.close_context_definition();

    // Create HField and ZField instances to read the data into
    HField mask(ModelArray::Type::H);
    mask.resize();
    VertexField coordinates(ModelArray::Type::VERTEX);
    coordinates.resize();
    DGField hice(ModelArray::Type::DG);
    hice.resize();

    // Setup ModelState with field above
    // FIXME: VertexField and DGField not read correctly
    ModelState state = { {
                             { maskName, mask },
                             { coordsName, coordinates },
                             { hiceName, hice },
                         },
        {} };

    // Check calendar step is zero initially
    REQUIRE(xiosHandler.getCalendarStep() == 0);

    // Check the input file exists
    REQUIRE(std::filesystem::exists(filename));

    // Deduce the local lengths of the two dimensions
    const size_t nx = ModelArray::definedDimensions.at(ModelArray::Dimension::X).localLength;
    const size_t ny = ModelArray::definedDimensions.at(ModelArray::Dimension::Y).localLength;

    // Simulate 4 iterations (timesteps)
    metadata.setTime(xiosHandler.getCalendarStart());
    for (int ts = 1; ts <= 4; ts++) {
        // Update the current timestep and verify it's updated in XIOS
        metadata.incrementTime(timestep);
        REQUIRE(xiosHandler.getCalendarStep() == ts);
        // Check that the fields contain the expected data
        ModelState state = grid.getModelState(filename, metadata);
        for (auto& entry : state.data) {
            if (entry.first == maskName) {
                for (size_t j = 0; j < ny; ++j) {
                    for (size_t i = 0; i < nx; ++i) {
                        REQUIRE(entry.second(i, j) == doctest::Approx(j >= 1 ? 1.0 : 0.0));
                    }
                }
            } else if (entry.first == coordsName) {
                for (size_t j = 0; j < ny + 1; ++j) {
                    for (size_t i = 0; i < nx + 1; ++i) {
                        REQUIRE(coordinates.components({ i, j })[0] == doctest::Approx(i));
                        REQUIRE(coordinates.components({ i, j })[1] == doctest::Approx(j));
                    }
                }
            } else if (entry.first == hiceName) {
                for (size_t j = 0; j < ny; ++j) {
                    for (size_t i = 0; i < nx; ++i) {
                        for (size_t d = 0; d < DG; ++d) {
                            float expected = 1.0 * (d + DG * (i + nx * j));
                            REQUIRE(hice.components({ i, j })[d] == doctest::Approx(expected));
                        }
                    }
                }
            }
        }
    }
    xiosHandler.context_finalize();
    Finalizer::finalize();
}

}
