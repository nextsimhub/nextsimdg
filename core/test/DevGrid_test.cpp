/*!
 * @file DevGrid_example.cpp
 *
 * @date Jan 7, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "include/DevGrid.hpp"
#include "include/DevGridIO.hpp"
#include "include/ElementData.hpp"
#include "include/IStructure.hpp"

#include <cstdio>
#include <fstream>

const std::string filename = "DevGrid_test.nc";

namespace Nextsim {

TEST_CASE("Write out a DevGrid restart file", "[DevGrid]")
{
    DevGrid grid;
    grid.init("");
    grid.setIO(new DevGridIO(grid));
    // Fill in the data. It is not real data.
    grid.resetCursor();
    int nx = DevGrid::nx;
    int ny = DevGrid::nx;
    double yFactor = 0.0001;
    double xFactor = 0.01;

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            if (grid.validCursor()) {
                double fractional = j * yFactor + i * xFactor;
                grid.cursorData() = PrognosticGenerator()
                                        .hice(1 + fractional)
                                        .cice(2 + fractional)
                                        .sst(3 + fractional)
                                        .sss(4 + fractional)
                                        .hsnow(5 + fractional)
                                        .tice({ -(1. + fractional) });
                grid.incrCursor();
            }
        }
    }

    grid.dump(filename);

    DevGrid grid2;
    grid2.init("");
    grid2.setIO(new DevGridIO(grid2));

    grid2.cursor = 0;

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            if (grid2.cursor) {
                *(grid2.cursor) = PrognosticGenerator();
                ++grid2.cursor;
            }
        }
    }

    grid2.cursor = 0;
    int targetIndex = 7 * DevGrid::nx + 3;
    for (int i = 0; i < targetIndex; ++i) {
        ++grid2.cursor;
    }
    if (!grid2.cursor) {
        FAIL("Invalid cursor value of " << targetIndex);
    }

    double unInitIce = grid2.cursor->iceThickness();
    REQUIRE(unInitIce == 0.);

    grid2.init(filename);

    grid2.cursor = 0;
    for (int i = 0; i < targetIndex; ++i) {
        ++grid2.cursor;
    }
    if (!grid2.cursor) {
        FAIL("Invalid cursor value of " << targetIndex);
    }

    REQUIRE(grid2.cursor->iceThickness() != 0);
    REQUIRE(grid2.cursor->iceThickness() > 1);
    REQUIRE(grid2.cursor->iceThickness() < 2);
    REQUIRE(grid2.cursor->iceThickness() == 1.0703);
    REQUIRE(grid2.cursor->iceThickness() != unInitIce);

    REQUIRE(grid2.cursor->iceTemperature(0) < -1);
    REQUIRE(grid2.cursor->iceTemperature(0) > -2);

    std::remove(filename.c_str());
}

const std::string stateFilename = "modelState.test.nc";

TEST_CASE("Write and read a ModelState-based DevGrid restart file", "[DevGrid]")
{
    DevGrid grid;
    grid.init("");
    grid.setIO(new DevGridIO(grid));

    // Fill in the data. It is not real data.
    size_t nx = DevGrid::nx;
    size_t ny = DevGrid::nx;
    double yFactor = 0.01;
    double xFactor = 0.0001;

    ModelArray::setDimensions(ModelArray::Type::H, {nx, ny});
    ModelArray::setDimensions(ModelArray::Type::Z, {nx, ny, 1});

    HField fractional(ModelArray::Type::H, "");
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            fractional(i, j) = j * yFactor + i * xFactor;
        }
    }

    HField hice = fractional + 1;
    HField cice = fractional + 2;
    HField sst = fractional + 3;
    HField sss = fractional + 4;
    HField hsnow = fractional + 5;

    HField ticeValue = -(fractional + 1);
    ZField tice = ModelArray::ZField("tice");
    tice.setData(ticeValue);

    ModelState state = {
            {"hice", hice},
            {"cice", cice},
            {"sst", sst},
            {"sss", sss},
            {"hsnow", hsnow},
            {"tice", tice},
    };

    grid.dumpModelState(state, stateFilename);

    ModelArray::setDimensions(ModelArray::Type::H, {1, 1});
    REQUIRE(ModelArray::dimensions(ModelArray::Type::H)[0] == 1);
    DevGrid gridIn;
    size_t targetX = 3;
    size_t targetY = 7;

    gridIn.setIO(new DevGridIO(grid));
    ModelState ms = gridIn.getModelState(stateFilename);

    REQUIRE(ModelArray::dimensions(ModelArray::Type::H)[0] == DevGrid::nx);
    REQUIRE(ModelArray::dimensions(ModelArray::Type::H)[1] == DevGrid::nx);
    REQUIRE(ms.at("hice")(targetX, targetY) != 0);
    REQUIRE(ms.at("hice")(targetX, targetY) > 1);
    REQUIRE(ms.at("hice")(targetX, targetY) < 2);
    REQUIRE(ms.at("hice")(targetX, targetY) == 1.0703);

    ZField ticeIn = ms.at("tice");

    REQUIRE(ticeIn.dimensions()[2] == 1);
    REQUIRE(ticeIn(targetX, targetY, 0) == -1.0703);

    std::remove(stateFilename.c_str());

}


}
