/*!
 * @file DevGrid_example.cpp
 *
 * @date Jan 7, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "include/DevGrid.hpp"
#include "include/ElementData.hpp"
#include "include/IStructure.hpp"
#include "include/ModuleLoader.hpp"

#include <cstdio>
#include <fstream>

const std::string filename = "DevGrid_test.nc";

namespace Nextsim {

TEST_CASE("Write out a DevGrid restart file", "[DevGrid]")
{
    ModuleLoader::getLoader().setAllDefaults();

    DevGrid grid;
    // Fill in the data. It is not real data.
    grid.resetCursor();
    int nx = DevGrid::nx;
    int ny = DevGrid::nx;
    double yFactor = 0.01;
    double xFactor = 0.0001;

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            if (grid.validCursor()) {
                double fractional = j * yFactor + i * xFactor;
                grid.cursorData() = PrognosticData::generate(1 + fractional, 2 + fractional,
                    3 + fractional, 4 + fractional, 5 + fractional, { 0., 0., 0. });
                grid.incrCursor();
            }
        }
    }

    grid.dump(filename);

    DevGrid grid2;

    grid2.resetCursor();

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            if (grid2.validCursor()) {
                grid2.cursorData() = PrognosticData::generate(0., 0, 0, 0, 0, { 0., 0., 0. });
                grid2.incrCursor();
            }
        }
    }

    grid2.resetCursor();
    int targetIndex = 7 * DevGrid::nx + 3;
    for (int i = 0; i < targetIndex; ++i) {
        grid2.incrCursor();
    }
    if (!grid2.validCursor()) {
        FAIL("Invalid cursor value of " << targetIndex);
    }

    double unInitIce = grid2.cursorData().iceThickness();
    REQUIRE(unInitIce == 0.);

    grid2.init(filename);

    grid2.resetCursor();
    for (int i = 0; i < targetIndex; ++i) {
        grid2.incrCursor();
    }
    if (!grid2.validCursor()) {
        FAIL("Invalid cursor value of " << targetIndex);
    }

    REQUIRE(grid2.cursorData().iceThickness() != 0);
    REQUIRE(grid2.cursorData().iceThickness() > 1);
    REQUIRE(grid2.cursorData().iceThickness() < 2);
    REQUIRE(grid2.cursorData().iceThickness() == 1.0703);
    REQUIRE(grid2.cursorData().iceThickness() != unInitIce);

    std::remove(filename.c_str());
}
}
