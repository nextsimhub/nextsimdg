/*!
 * @file RectGrid_test.cpp
 *
 * @date Feb 8, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "include/RectangularGrid.hpp"
#include "include/RectGridIO.hpp"
#include "include/ElementData.hpp"
#include "include/IStructure.hpp"

#include <cstdio>
#include <fstream>

const std::string filename = "RectGrid_test.nc";

namespace Nextsim {

TEST_CASE("Write out a RectangularGrid restart file", "[RectangularGrid]")
{
    int nx = 25;
    int ny = 35;
    int nLayers = 2;

    RectangularGrid::GridDimensions dims = { nx, ny, nLayers };

    RectangularGrid grid(dims);
    grid.setIO(new RectGridIO(grid));
    // Fill in the data. It is not real data.
    grid.resetCursor();

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
                                        .tice({ -(1. + fractional), -(2. +fractional) });
                grid.incrCursor();
            }
        }
    }

    grid.dump(filename);

    RectangularGrid grid2;
    grid2.setIO(new RectGridIO(grid));

    int targetIndex = 17 * ny + 33;
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
    REQUIRE(grid2.cursor->iceThickness() == 1.1733);

    REQUIRE(grid2.cursor->iceTemperature(1) < -2);
    REQUIRE(grid2.cursor->iceTemperature(1) > -3);
    REQUIRE(grid2.cursor->iceTemperature(1) == Approx(-2.1733).epsilon(1e-14));

    std::remove(filename.c_str());
}
} /* namespace Nextsim */
