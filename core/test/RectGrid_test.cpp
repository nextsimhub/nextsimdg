/*!
 * @file RectGrid_test.cpp
 *
 * @date 24 Sep 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifdef USE_MPI
#include <doctest/extensions/doctest_mpi.h>
#else
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#endif

#include "include/CommonRestartMetadata.hpp"
#include "include/IStructure.hpp"
#include "include/NZLevels.hpp"
#include "include/RectGridIO.hpp"
#include "include/RectangularGrid.hpp"
#include "include/gridNames.hpp"

#include <cstdio>
#include <fstream>

const std::string test_files_dir = TEST_FILES_DIR;
const std::string filename = test_files_dir + "/RectGrid_test.nc";
#ifdef USE_MPI
const std::string filename_parallel = test_files_dir + "/RectGrid_test_parallel.nc";
const std::string partition_filename = test_files_dir + "/partition_metadata_3.nc";
#endif
const std::string date_string = "2000-01-01T00:00:00Z";

namespace Nextsim {
TEST_SUITE_BEGIN("RectGrid");
#ifdef USE_MPI
// Number of ranks should not be hardcoded here
MPI_TEST_CASE("Write and read a ModelState-based RectGrid restart file", 3)
#else
TEST_CASE("Write and read a ModelState-based RectGrid restart file")
#endif
{
    RectangularGrid grid;
    grid.setIO(new RectGridIO(grid));

    // Fill in the data. It is not real data.
    size_t nx = 5;
    size_t ny = 7;
    double yFactor = 0.01;
    double xFactor = 0.0001;

    // Create data for reference file
    NZLevels::set(1);
    ModelArray::setDimensions(ModelArray::Type::H, { nx, ny });
    ModelArray::setDimensions(ModelArray::Type::Z, { nx, ny, NZLevels::get() });

    HField fractional(ModelArray::Type::H);
    HField mask(ModelArray::Type::H);
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            fractional(i, j) = j * yFactor + i * xFactor;
            mask(i, j)
                = (i - nx / 2) * (i - nx / 2) + (j - ny / 2) * (j - ny / 2) > (nx * ny) ? 0 : 1;
        }
    }

    HField hice = fractional + 1;
    HField cice = fractional + 2;
    HField sst = fractional + 3;
    HField sss = fractional + 4;
    HField hsnow = fractional + 5;

    HField ticeValue = -(fractional + 1);
    ZField tice = ModelArray::ZField();
    tice.setData(ticeValue);

    ModelState state = { {
                             { "mask", mask },
                             { "hice", hice },
                             { "cice", cice },
                             { "hsnow", hsnow },
                             { "tice", tice },
                         },
        {} };

    ModelMetadata metadata;
    metadata.setTime(TimePoint(date_string));
    // Use x & y coordinates
    ModelArray x(ModelArray::Type::H);
    ModelArray y(ModelArray::Type::H);
    // Use an anisotropic grid so we can differentiate the dimensions.
    double dx = 25.;
    double dy = 35.;
    for (int j = 0; j < ny; ++j) {
        double yy = j * dy;
        for (int i = 0; i < nx; ++i) {
            double xx = i * dx;
            x(i, j) = xx;
            y(i, j) = yy;
        }
    }
    // Use a temporary state to set the coordinates
    ModelState coordState = { {
                                  { xName, x },
                                  { yName, y },
                              },
        {} };
    metadata.extractCoordinates(coordState);
    // Then immediately extract them to the output state
    metadata.affixCoordinates(state);

// Write reference file
#ifdef USE_MPI
    // Create subcommunicator with only first rank
    metadata.setMpiMetadata(test_comm);
    int colour = MPI_UNDEFINED, key = 0;
    MPI_Comm rank0Comm;

    if (metadata.mpiMyRank == 0) {
        colour = 0;
    }
    MPI_Comm_split(test_comm, colour, key, &rank0Comm);

    // Write reference file serially on first MPI rank
    if (metadata.mpiMyRank == 0) {
        metadata.setMpiMetadata(rank0Comm);
        metadata.globalExtentX = nx;
        metadata.globalExtentY = ny;
        metadata.localCornerX = 0;
        metadata.localCornerY = 0;
        metadata.localExtentX = nx;
        metadata.localExtentY = ny;
        grid.dumpModelState(state, metadata, filename);
        MPI_Comm_free(&rank0Comm);
    }
#else
    grid.dumpModelState(state, metadata, filename);
#endif

    // Reset dimensions so it is possible to check if they
    // are read correctly from refeence file
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });
    ModelArray::setDimensions(ModelArray::Type::Z, { 1, 1, 1 });
    REQUIRE(ModelArray::dimensions(ModelArray::Type::H)[0] == 1);
    RectangularGrid gridIn;
    size_t targetX = 1;
    size_t targetY = 2;

    // Read reference file
    gridIn.setIO(new RectGridIO(grid));
#ifdef USE_MPI
    ModelMetadata metadataIn(partition_filename, test_comm);
    metadataIn.setTime(TimePoint(date_string));
    ModelState ms = gridIn.getModelState(filename, metadataIn);
#else
    ModelState ms = gridIn.getModelState(filename);
#endif

// Check if dimensions and data from reference file are correct
#ifdef USE_MPI
    REQUIRE(ModelArray::dimensions(ModelArray::Type::H)[0] == metadataIn.localExtentX);
    REQUIRE(ModelArray::dimensions(ModelArray::Type::H)[1] == metadataIn.localExtentY);
    REQUIRE(ModelArray::dimensions(ModelArray::Type::Z)[0] == metadataIn.localExtentX);
    REQUIRE(ModelArray::dimensions(ModelArray::Type::Z)[1] == metadataIn.localExtentY);
#else
    REQUIRE(ModelArray::dimensions(ModelArray::Type::H)[0] == nx);
    REQUIRE(ModelArray::dimensions(ModelArray::Type::H)[1] == ny);
    REQUIRE(ModelArray::dimensions(ModelArray::Type::Z)[0] == nx);
    REQUIRE(ModelArray::dimensions(ModelArray::Type::Z)[1] == ny);
#endif
#ifdef USE_MPI
    REQUIRE(ms.data.at("hice")(targetX, targetY)
        == 1.0201 + metadataIn.localCornerY * yFactor + metadataIn.localCornerX * xFactor);
#else
    REQUIRE(ms.data.at("hice")(targetX, targetY) == 1.0201);
#endif

    ZField ticeIn = ms.data.at("tice");

    REQUIRE(ticeIn.dimensions()[2] == 1);
#ifdef USE_MPI
    REQUIRE(ticeIn(targetX, targetY, 0U)
        == -1.0201 - metadataIn.localCornerY * yFactor - metadataIn.localCornerX * xFactor);
#else
    REQUIRE(ticeIn(targetX, targetY, 0U) == -1.0201);
#endif

    // Check that the coordinates have been correctly written and read
    REQUIRE(ms.data.count(xName) > 0);
    REQUIRE(ms.data.count(yName) > 0);
#ifdef USE_MPI
    REQUIRE(ms.data.at(xName)(1, 0) == dx * (metadataIn.localCornerX + 1));
    REQUIRE(ms.data.at(xName)(0, 1) == dx * metadataIn.localCornerX);
    REQUIRE(ms.data.at(yName)(0, 1) == dy * (metadataIn.localCornerY + 1));
#else
    REQUIRE(ms.data.at(xName)(1, 0) == dx);
    REQUIRE(ms.data.at(xName)(0, 1) == 0);
    REQUIRE(ms.data.at(yName)(0, 1) == dy);
#endif

// Write file in parallel so it can be compared with one written serially
#ifdef USE_MPI
    gridIn.dumpModelState(ms, metadataIn, filename_parallel);
#endif

    std::remove(filename.c_str());
}
TEST_SUITE_END();

} /* namespace Nextsim */
