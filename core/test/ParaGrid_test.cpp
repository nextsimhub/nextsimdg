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

#include <fstream>

const std::string filename = "ParaGrid_test.nc";

static const int DG = 3;
static const int DGSTRESS = 6;
static const int CG = 2;

namespace Nextsim {

TEST_CASE("Write and read a ModelState-based ParaGrid restart file", "[ParametricGrid]")
{
    ParametricGrid grid;
    ParaGridIO pio(grid);

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
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            fractional(i, j) = j * yFactor + i * xFactor;
            mask(i, j) = (i - nx / 2)*(i - nx/2) + (j - ny / 2) * (j - ny / 2)  > (nx * ny) ? 0 : 1;
            for (size_t d = 0; d < DG; ++d) {
                fractionalDG.components({i, j})[d] = fractional(i, j) + d;
            }
        }
    }
}

}
