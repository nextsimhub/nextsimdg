/*!
 * @file ParaGrid_test.cpp
 *
 * @date Oct 27, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "include/CommonRestartMetadata.hpp"
#include "include/ParametricGrid.hpp"
#include "include/ParaGridIO.hpp"
#include "include/gridNames.hpp"

#include <cstdio>
#include <fstream>

const std::string filename = "paraGrid_test.nc";

static const int DG = 3;
static const int DGSTRESS = 6;
static const int CG = 2;

namespace Nextsim {

TEST_CASE("Write and read a ModelState-based ParaGrid restart file", "[ParametricGrid]")
{
    std::FILE* lun = std::fopen(filename.c_str(), "r");
    if (lun != NULL) {
        std::remove(filename.c_str());
    }
    std::fclose(lun);

    ParametricGrid grid;
    ParaGridIO writeIO(grid);
    grid.setIO(&writeIO);

    // Set the dimension lengths
    size_t nx = 25;
    size_t ny = 15;
    size_t nz = 3;
    size_t nxcg = CG * nx + 1;
    size_t nycg = CG * ny + 1;

    double yFactor = 0.01;
    double xFactor = 0.0001;

    ModelArray::setDimension(ModelArray::Dimension::X, nx);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny);
    ModelArray::setDimension(ModelArray::Dimension::Z, nz);
    ModelArray::setDimension(ModelArray::Dimension::XCG, nxcg);
    ModelArray::setDimension(ModelArray::Dimension::YCG, nycg);

    ModelArray::setNComponents(ModelArray::Type::DG, DG);
    ModelArray::setNComponents(ModelArray::Type::DGSTRESS, DGSTRESS);

    HField fractional(ModelArray::Type::H);
    DGField fractionalDG(ModelArray::Type::DG);
    HField mask(ModelArray::Type::H);
    fractional.resize();
    fractionalDG.resize();
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            fractional(i, j) = j * yFactor + i * xFactor;
            mask(i, j) = (i - nx / 2)*(i - nx/2) + (j - ny / 2) * (j - ny / 2)  > (nx * ny) ? 0 : 1;
            for (size_t d = 0; d < DG; ++d) {
                fractionalDG.components({i, j})[d] = fractional(i, j) + d;
            }
        }
    }

    DGField hice = fractionalDG + 10;
    DGField cice = fractionalDG + 20;
    DGField hsnow = fractionalDG + 30;
    ZField tice(ModelArray::Type::Z);
    tice.resize();
    for (size_t i = 0; i < ModelArray::size(ModelArray::Type::H); ++i) {
        for (size_t k = 0; k < nz; ++k) {
            tice.zIndexAndLayer(i, k) = fractional[i] + 40 + k;
        }
    }

    ModelState state = {{
            { maskName, mask },
            { hiceName, hice },
            { ciceName, cice },
            { hsnowName, hsnow },
            { ticeName, tice },
    }, {}};

    ModelMetadata metadata;
    metadata.setTime(TimePoint("2000-01-01T00:00:00Z"));

    grid.dumpModelState(state, metadata, filename, true);

    lun = std::fopen(filename.c_str(), "r");
    REQUIRE(lun != NULL);
    REQUIRE(std::fclose(lun) == 0);

    // Reset the array dimensions to make sure that the read function gets them correct
    ModelArray::setDimension(ModelArray::Dimension::X, 1);
    ModelArray::setDimension(ModelArray::Dimension::Y, 1);
    ModelArray::setDimension(ModelArray::Dimension::Z, 1);
    ModelArray::setDimension(ModelArray::Dimension::XCG, 1);
    ModelArray::setDimension(ModelArray::Dimension::YCG, 1);
    // In the full model numbers of DG components are set at compile time, so they are not reset
    REQUIRE(ModelArray::nComponents(ModelArray::Type::DG) == DG);
    
    ParametricGrid gridIn;
    ParaGridIO readIO(gridIn);
    gridIn.setIO(&readIO);

    ModelState ms = gridIn.getModelState(filename);

    REQUIRE(ModelArray::dimensions(ModelArray::Type::Z)[0] == nx);
    REQUIRE(ModelArray::dimensions(ModelArray::Type::Z)[1] == ny);
    REQUIRE(ModelArray::dimensions(ModelArray::Type::Z)[2] == nz);

    REQUIRE( ms.data.size() == state.data.size() );

    ModelArray& ticeRef = ms.data.at(ticeName);
    REQUIRE(ModelArray::nDimensions(ModelArray::Type::Z) == 3);
    REQUIRE(ticeRef.getType() == ModelArray::Type::Z);
    REQUIRE(ticeRef.nDimensions() == 3);
    REQUIRE(ticeRef.dimensions()[0] == nx);
    REQUIRE(ticeRef.dimensions()[1] == ny);
    REQUIRE(ticeRef.dimensions()[2] == nz);

    ModelArray& hiceRef = ms.data.at(hiceName);
    REQUIRE(hiceRef.nDimensions() == 2);
    REQUIRE(hiceRef.dimensions()[0] == nx);
    REQUIRE(hiceRef.dimensions()[1] == ny);
    REQUIRE(ModelArray::nComponents(ModelArray::Type::DG) == DG);
    REQUIRE(hiceRef.nComponents() == DG);

    REQUIRE(ticeRef(12, 14, 1) == tice(12, 14, 1));

    std::remove(filename.c_str());
}
}
