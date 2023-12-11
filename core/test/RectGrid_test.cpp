/*!
 * @file RectGrid_test.cpp
 *
 * @date Feb 8, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifdef USE_MPI
#include <doctest/extensions/doctest_mpi.h>
#else
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#endif


#include "include/CommonRestartMetadata.hpp"
#include "include/NZLevels.hpp"
#include "include/RectangularGrid.hpp"
#include "include/RectGridIO.hpp"
#include "include/IStructure.hpp"

#include <cstdio>
#include <fstream>

const std::string filename = "RectGrid_test.nc";
const std::string partition_filename = "partition_metadata_1.nc";

namespace Nextsim {
TEST_SUITE_BEGIN("RectGrid");
#ifdef USE_MPI
// Number of ranks should not be hardcoded here
MPI_TEST_CASE("Write and read a ModelState-based RectGrid restart file", 1)
#else
TEST_CASE("Write and read a ModelState-based RectGrid restart file")
#endif
{
    RectangularGrid grid;
    grid.setIO(new RectGridIO(grid));

    // Fill in the data. It is not real data.
    size_t nx = 25;
    size_t ny = 15;
    double yFactor = 0.01;
    double xFactor = 0.0001;

    NZLevels::set(1);
    ModelArray::setDimensions(ModelArray::Type::H, { nx, ny });
    ModelArray::setDimensions(ModelArray::Type::Z, { nx, ny, NZLevels::get() });

    HField fractional(ModelArray::Type::H);
    HField mask(ModelArray::Type::H);
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            fractional(i, j) = j * yFactor + i * xFactor;
            mask(i, j) = (i - nx / 2)*(i - nx/2) + (j - ny / 2) * (j - ny / 2)  > (nx * ny) ? 0 : 1;
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

    ModelState state = {{
        { "mask", mask },
        { "hice", hice },
        { "cice", cice },
        { "sst", sst },
        { "sss", sss },
        { "hsnow", hsnow },
        { "tice", tice },
    }, {}};

    ModelMetadata metadata;
    metadata.setTime(TimePoint("2000-01-01T00:00:00Z"));

#ifdef USE_MPI
    metadata.setMpiMetadata(test_comm);
    metadata.globalExtentX = nx;
    metadata.globalExtentY = ny;
    metadata.localCornerX = 0;
    metadata.localCornerY = 0;
    metadata.localExtentX = nx;
    metadata.localExtentY = ny;
#endif
    grid.dumpModelState(state, metadata, filename);

    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });
    REQUIRE(ModelArray::dimensions(ModelArray::Type::H)[0] == 1);
    RectangularGrid gridIn;
    size_t targetX = 3;
    size_t targetY = 7;

    gridIn.setIO(new RectGridIO(grid));
#ifdef USE_MPI
    ModelState ms = gridIn.getModelState(filename, partition_filename, metadata);
#else
    ModelState ms = gridIn.getModelState(filename);
#endif

    REQUIRE(ModelArray::dimensions(ModelArray::Type::H)[0] == nx);
    REQUIRE(ModelArray::dimensions(ModelArray::Type::H)[1] == ny);
    REQUIRE(ms.data.at("hice")(targetX, targetY) != 0);
    REQUIRE(ms.data.at("hice")(targetX, targetY) > 1);
    REQUIRE(ms.data.at("hice")(targetX, targetY) < 2);
    REQUIRE(ms.data.at("hice")(targetX, targetY) == 1.0703);

    ZField ticeIn = ms.data.at("tice");

    REQUIRE(ticeIn.dimensions()[2] == 1);
    REQUIRE(ticeIn(targetX, targetY, 0U) == -1.0703);

    std::remove(filename.c_str());
}
TEST_SUITE_END();

} /* namespace Nextsim */
