#include "ParametricMesh.hpp"

#include <fstream>
#include <iostream>

namespace Nextsim {

void ParametricMesh::readmesh(std::string fname)
{
    reset();

    std::ifstream IN(fname.c_str());
    if (IN.fail()) {
        std::cerr << "ParametricMesh :: Could not open mesh file " << fname << std::endl;
        abort();
    }

    std::string status;
    IN >> status;

    if (status != "ParametricMesh") {
        std::cerr << "ParametricMesh :: Wrong file format " << fname << "\t" << status << std::endl
                  << "'ParametricMesh' expected" << std::endl;
        abort();
    }

    std::string version;
    IN >> version;

    if (statuslog > 0)
        std::cout << "ParametricMesh :: Reading " << fname << " V" << version << std::endl;

    if ((version != "1.0") && (version != "2.0")) {
        std::cerr << "ParametricMesh :: Wrong file format version " << fname << std::endl;
        abort();
    }

    IN >> nx >> ny;

    if (statuslog > 0)
        std::cout << "ParametricMesh :: Reading mesh with " << nx << " * " << ny << " elements"
                  << std::endl;

    if ((nx < 1) || (ny < 1)) {
        std::cerr << "ParametricMesh :: Wrong mesh dimensions (nx,ny) << " << nx << " " << ny
                  << std::endl;
        abort();
    }

    // set number of elements & nodes
    nelements = nx * ny;
    nnodes = (nx + 1) * (ny + 1);
    vertices.resize(nnodes, 2);

    for (size_t i = 0; i < nnodes; ++i) {
        IN >> vertices(i, 0) >> vertices(i, 1);
        if (IN.eof()) {
            std::cerr << "ParametricMesh :: Unexpected eof << " << fname << std::endl;
            abort();
        }
    }

    // Boundary
    if (version == "1.0") // all four boundaries are dirichlet
    {
        for (size_t i = 0; i < nx; ++i) // lower
            dirichlet[0].push_back(i);

        for (size_t i = 0; i < ny; ++i) // right
            dirichlet[1].push_back(i * nx + nx - 1);

        for (size_t i = 0; i < nx; ++i) // upper
            dirichlet[2].push_back(i + nx * (ny - 1));

        for (size_t i = 0; i < ny; ++i) // left
            dirichlet[3].push_back(i * nx);

        landmask.resize(nx * ny, true); // set landmask
    } else if (version == "2.0") // landmask, dirichlet and periodic boundary as additional lists
    {
        IN >> status;

        if (status != "landmask") {
            std::cerr << "V2.0 Expecting landmask information after nodes." << std::endl
                      << "\tlandmask ne" << std::endl
                      << "where ne is the number of elements. Should match nx * ny" << std::endl
                      << "If all is land, just provide " << std::endl
                      << "   landmask 0" << std::endl;
            abort();
        }
        size_t ne;
        IN >> ne;

        if (ne == 0) {
            landmask.resize(nx * ny, true);
        } else {
            assert(ne == nx * ny);
            landmask.resize(nx * ny, false);
            bool lm;
            for (size_t i = 0; i < nx * ny; ++i) {
                IN >> lm;

                if (lm)
                    landmask[i] = true;

                if (IN.eof()) {
                    std::cerr << "ParametricMesh :: Unexpected eof << " << fname << std::endl;
                    abort();
                }
            }
        }

        IN >> status;

        if (status != "dirichlet") {
            std::cerr << "V2.0 Expecting Dirichlet information after list of elements" << std::endl
                      << "\tdirichlet nd" << std::endl
                      << "where nd is the number of dirichlet edges" << std::endl
                      << "Then, for each edge we expect a line with two entries: id-of-element "
                         "[0,1,2,3] for the edge"
                      << std::endl;
            abort();
        }
        size_t nd;
        IN >> nd;

        if (statuslog > 0)
            std::cout << "reading " << nd << " dirichlet segments" << std::endl;

        for (size_t i = 0; i < nd; ++i) {
            size_t n0, n1;
            IN >> n0 >> n1; // read the element and the side
            assert(n0 < nx * ny);
            assert(n0 >= 0);

            dirichlet[n1].push_back(n0);

            if (IN.eof()) {
                std::cerr << "ParametricMesh :: Unexpected eof << " << fname << std::endl;
                abort();
            }
        }

        IN >> status;
        if (status != "periodic") {
            std::cerr << "Expecting Periodic information after Dirichlet" << std::endl
                      << "\tperiodic nd" << std::endl
                      << "where nd is the number of periodic segments" << std::endl;
            std::cerr << "Instead, got: " << status << std::endl;
            abort();
        }
        IN >> nd;
        periodic.clear();
        periodic.resize(nd);
        for (size_t i = 0; i < nd; ++i) {
            periodic[i].clear();
            size_t nind;
            IN >> nind;
            if (statuslog > 0)
                std::cout << "reading periodic segment " << i << " with " << nind << " entries"
                          << std::endl;

            for (size_t j = 0; j < nind; ++j) {
                size_t n0, n1, n2;
                IN >> n0 >> n1 >> n2; // read the elements and the side
                assert(n2 == 0 || n2 == 1);
                if (n2 == 0) // X-edge, bottom / top
                {
                    periodic[i].push_back(std::array<size_t, 4>({ n2, n0, n1, n0 }));
                } else if (n2 == 1) // Y-edge, left / right
                {
                    size_t ix = n0 % nx;
                    size_t iy = n0 / nx;
                    periodic[i].push_back(
                        std::array<size_t, 4>({ n2, n0, n1, iy * (nx + 1) + ix }));
                } else
                    abort();
            }
        }
    }

    IN.close();

    if (statuslog > 0) {
        std::cout << "ParametricMesh :: read mesh file " << fname << std::endl
                  << "             nx,ny = " << nx << " , " << ny << std::endl
                  << "             " << nelements << " elements,  " << nnodes << " nodes"
                  << std::endl;
    }
}

void ParametricMesh::coordinatesFromModelArray(const ModelArray& coord1, const ModelArray& coord2)
{
    // Fill in the array sizes from the ModelArray dimensions
    nx = ModelArray::size(ModelArray::Dimension::X);
    ny = ModelArray::size(ModelArray::Dimension::Y);
    nelements = nx * ny;
    nnodes = (nx + 1) * (ny + 1);
    vertices.resize(nnodes, 2);
    for (size_t idx = 0; idx < nnodes; ++idx) {
        vertices(idx, 0) = coord1(idx);
        vertices(idx, 1) = coord2(idx);
    }
}

void ParametricMesh::dirichletFromMask(const ModelArray& mask)
{
    // BOTTOM, RIGHT, TOP, LEFT
    const std::array<size_t, N_EDGE> startX = {0, 0, 0, 1};
    const std::array<size_t, N_EDGE> stopX = {nx, nx - 1, nx, nx};
    const std::array<size_t, N_EDGE> startY = {1, 0, 0, 0};
    const std::array<size_t, N_EDGE> stopY = {ny, ny, ny - 1, ny};
    const std::array<int, N_EDGE> deltaX = {0, 1, 0, -1};
    const std::array<int, N_EDGE> deltaY = {-1, 0, 1, 0};

    // Loop over edges
    for (auto edge = edges.begin(); edge != edges.end(); ++edge) {
        for (size_t j = startY[*edge]; j < stopY[*edge]; ++j) {
            for (size_t i = startX[*edge]; i < stopX[*edge]; ++i) {
                if (!mask(i, j)) continue;
                // mask(i, j) is ocean. Check the appropriate neighbour
                if (!mask(i + deltaX[*edge], j + deltaY[*edge])) {
                    dirichlet[*edge].push_back(mask.indexFromLocation({i, j}));
                }
            }
        }
    }
}

void ParametricMesh::dirichletFromEdge(const ModelArray& mask, Edge edge)
{
    // BOTTOM, RIGHT, TOP, LEFT
    const std::array<size_t, N_EDGE> startX = {0, nx, 0, 0};
    const std::array<size_t, N_EDGE> stopX = {nx, nx, nx, 0};
    const std::array<size_t, N_EDGE> startY = {0, 0, ny, 0};
    const std::array<size_t, N_EDGE> stopY = {0, ny, ny, ny};

    for (size_t j = startY[edge]; j < stopY[edge]; ++j) {
        for (size_t i = startX[edge]; i < stopX[edge]; ++i) {
            if (mask(i, j)) {
                dirichlet[edge].push_back(mask.indexFromLocation({i, j}));
            }
        }
    }
}

void ParametricMesh::sortDirichlet()
{
    for (auto edge = edges.begin(); edge != edges.end(); ++edge) {
        std::sort(dirichlet[*edge].begin(), dirichlet[*edge].end());
    }
}
/*!
 * returns minimum mesh size.
 *
 */
double ParametricMesh::hmin() const
{
    double hmin = 1.e99;
    for (size_t i = 0; i < nelements; ++i)
        hmin = std::min(hmin, h(i));
    return hmin;
}

/*!
 * returns are of domain
 */
double ParametricMesh::area() const
{
    double a = 0;
    for (size_t i = 0; i < nelements; ++i)
        a += area(i);
    return a;
}

}
