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

    // Create the ModelState to be written
    // Set the size of the arrays
    size_t snx = nx;
    size_t sny = ny;
    size_t snLayers = nLayers;
    ModelArray::Dimensions hdim = {snx, sny};
    ModelArray::setDimensions(ModelArray::Type::H, hdim);
    ModelArray::setDimensions(ModelArray::Type::U, hdim);
    ModelArray::setDimensions(ModelArray::Type::V, hdim);
    ModelArray::Dimensions zdim = {snx, sny, snLayers};
    ModelArray::setDimensions(ModelArray::Type::Z, zdim);

    // Create the data arrays, storing them in the ModelState
    ModelState ms;
    for (const std::string fieldName : { "hice", "cice", "hsnow", "sst", "sss"}) {
        ms[fieldName] = ModelArray::HField(fieldName);
    }
    std::string tice = "tice";
    ms[tice] = ModelArray::ZField(tice);

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            double fractional = j * yFactor + i * xFactor;
            ms["hice"](i, j) = (1 + fractional);
            ms["cice"](i, j) = (2 + fractional);
            ms["sst"](i, j) = (3 + fractional);
            ms["sss"](i, j) = (4 + fractional);
            ms["hsnow"](i, j) = (5 + fractional);
            ms[tice](i, j, 0) = (-(1. + fractional));
            ms[tice](i, j, 1) = (-(2. + fractional));
        }
    }

    RectangularGrid grid3;
    grid3.setIO(new RectGridIO(grid));
    grid3.dumpModelState(ms, filename);


// Read the file again, this time using the ModelState system
    RectangularGrid grid4;
    grid4.setIO(new RectGridIO(grid));
    ModelState state = grid4.getModelState(filename);
    ModelArray::Dimensions loc = {17, 33};

    std::string hice = "hice";
    REQUIRE(state[hice][loc] != 0);
    REQUIRE(state[hice][loc] > 1);
    REQUIRE(state[hice][loc] < 2);
    REQUIRE(state[hice][loc] == 1.1733);


    ModelArray::Dimensions zoc = {loc[0], loc[1], 1};
    REQUIRE(state[tice][zoc] < -2);
    REQUIRE(state[tice][zoc] > -3);
    REQUIRE(state[tice][zoc] == Approx(-2.1733).epsilon(1e-14));

    std::remove(filename.c_str());
}
} /* namespace Nextsim */
