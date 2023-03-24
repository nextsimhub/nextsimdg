/*!
 * @file DevGrid_example.cpp
 *
 * @date Jan 7, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "include/DevGrid.hpp"
#include "include/DevGridIO.hpp"
#include "include/IStructure.hpp"

#include <fstream>

const std::string filename = "example.nc";
const std::string stateFilename = "modelState.nc";

namespace Nextsim {

TEST_CASE("Write out a ModelState-based DevGrid restart file")
{
    DevGrid grid;
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
    WARN("File written to " << stateFilename);
}

TEST_CASE("Read the ModelState example file, if it exists")
{
    ModelArray::setDimensions(ModelArray::Type::H, {1, 1});
    REQUIRE(ModelArray::dimensions(ModelArray::Type::H)[0] == 1);
    DevGrid grid;
    size_t targetX = 3;
    size_t targetY = 7;

    grid.setIO(new DevGridIO(grid));
    ModelState ms = grid.getModelState(stateFilename);

    REQUIRE(ModelArray::dimensions(ModelArray::Type::H)[0] == DevGrid::nx);
    REQUIRE(ModelArray::dimensions(ModelArray::Type::H)[1] == DevGrid::nx);
    REQUIRE(ms.at("hice")(targetX, targetY) != 0);
    REQUIRE(ms.at("hice")(targetX, targetY) > 1);
    REQUIRE(ms.at("hice")(targetX, targetY) < 2);
    REQUIRE(ms.at("hice")(targetX, targetY) == 1.0703);
}
}
