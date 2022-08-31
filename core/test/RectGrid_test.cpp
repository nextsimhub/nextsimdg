/*!
 * @file RectGrid_test.cpp
 *
 * @date Feb 8, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "include/CommonRestartMetadata.hpp"
#include "include/RectangularGrid.hpp"
#include "include/RectGridIO.hpp"
#include "include/IStructure.hpp"

#include <cstdio>
#include <fstream>

const std::string filename = "RectGrid_test.nc";

namespace Nextsim {
TEST_CASE("Write and read a ModelState-based RectGrid restart file", "[DevGrid]")
{
    RectangularGrid grid;
    grid.setIO(new RectGridIO(grid));

    // Fill in the data. It is not real data.
    size_t nx = 25;
    size_t ny = 15;
    double yFactor = 0.01;
    double xFactor = 0.0001;

    ModelArray::setDimensions(ModelArray::Type::H, { nx, ny });
    ModelArray::setDimensions(ModelArray::Type::Z, { nx, ny, 1 });

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

    grid.dumpModelState(state, metadata, filename);

    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });
    REQUIRE(ModelArray::dimensions(ModelArray::Type::H)[0] == 1);
    RectangularGrid gridIn;
    size_t targetX = 3;
    size_t targetY = 7;

    gridIn.setIO(new RectGridIO(grid));
    ModelState ms = gridIn.getModelState(filename);

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


} /* namespace Nextsim */
