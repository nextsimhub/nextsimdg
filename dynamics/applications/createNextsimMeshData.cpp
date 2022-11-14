/*!
 * @file createNextsimMeshData.cpp
 *
 * @date Nov 3, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ParametricGrid.hpp"
#include "include/ParaGridIO.hpp"
#include "include/gridNames.hpp"

#include <cassert>
#include <string>

#define DGadvect 1
#define DGstress 3
#define CG 1

// Recreates the grid and initial data of benchmark_mehlmann_sasip_mevp.cpp as
// a NextsimDG netCDF restart file.

namespace Nextsim {
int makeMEVPFile(std::string filename)
{
    ParametricGrid grid;
    grid.setIO(new ParaGridIO(grid));

    // Set the dimension lengths
    size_t nx = 128;
    size_t ny = 128;
    size_t nz = 3;
    size_t nxcg = CG * nx + 1;
    size_t nycg = CG * ny + 1;

    ModelArray::setDimension(ModelArray::Dimension::X, nx);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny);
    ModelArray::setDimension(ModelArray::Dimension::Z, nz);
    ModelArray::setDimension(ModelArray::Dimension::XVERTEX, nx + 1);
    ModelArray::setDimension(ModelArray::Dimension::YVERTEX, ny + 1);
    ModelArray::setDimension(ModelArray::Dimension::XCG, nxcg);
    ModelArray::setDimension(ModelArray::Dimension::YCG, nycg);

    ModelArray::setNComponents(ModelArray::Type::DG, DGadvect);
    ModelArray::setNComponents(ModelArray::Type::DGSTRESS, DGstress);
    ModelArray::setNComponents(ModelArray::Type::VERTEX, ModelArray::nCoords);

    // Create the rectangle mesh, taken from createdistortedrectanglemesh.py
    VertexField coordinates(ModelArray::Type::VERTEX);
    coordinates.resize();
    double lx = 512000.;
    double ly = 512000.;
    for (size_t i = 0; i < ModelArray::definedDimensions.at(ModelArray::Dimension::XVERTEX).length; ++i) {
        for (size_t j = 0; j < ModelArray::definedDimensions.at(ModelArray::Dimension::YVERTEX).length; ++j) {
            // Swapped to match the results of the mesh file.
            coordinates.components({i, j})[1] = i * lx / nx;
            coordinates.components({i, j})[0] = j * ly / ny;
        }
    }

    DGField cice(ModelArray::Type::DG);
    DGField hice(ModelArray::Type::DG);
    cice.resize();
    hice.resize();
    double hiceMean = 0.3;
    double hiceAmpl = 0.005;
    double kx = 6e-5;
    double ky = 3e-5;
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            cice(i, j) = 1.0;
            hice(i, j) = hiceMean + 0.005 * (std::sin(kx * coordinates(i, j, 0UL)) + std::sin(ky * coordinates(i, j, 1UL)));
        }
    }

    ModelState state = {{
            { hiceName, hice },
            { ciceName, cice },
            { coordsName, coordinates },
    }, {}};

    ModelMetadata metadata;
    metadata.setTime(TimePoint("2000-01-01T00:00:00Z"));

    grid.dumpModelState(state, metadata, filename, true);

    return 0;
}
}

int main()
{
    std::string filename = "nextsimMEVPstart.nc";

    return Nextsim::makeMEVPFile(filename);
}
