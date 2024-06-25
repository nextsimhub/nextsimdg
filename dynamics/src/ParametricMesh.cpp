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
        std::cerr << "ParametricMesh :: Wrong file format " << fname << "\t" << status
                  << std::endl
                  << "'ParametricMesh' expected" << std::endl;
        abort();
    }

    std::string version;
    IN >> version;
    
    if (statuslog>0)
      std::cout << "ParametricMesh :: Reading " << fname << " V" << version  << std::endl;

    
    if ((version != "1.0") && (version != "2.0")) {
        std::cerr << "ParametricMesh :: Wrong file format version " << fname << std::endl;
        abort();
    }

    IN >> nx >> ny;

    if (statuslog>0)
      std::cout << "ParametricMesh :: Reading mesh with " << nx << " * " << ny << " elements" << std::endl;

    if ((nx < 1) || (ny < 1)) {
        std::cerr << "ParametricMesh :: Wrong mesh dimensions (nx,ny) << " << nx << " " << ny << std::endl;
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

      landmask.resize(nx*ny,true); // set landmask
    } else if (version == "2.0") // landmask, dirichlet and periodic boundary as additional lists
    {
        IN >> status;
        if (status != "landmask") {
	  std::cerr << "V2.0 Expecting landmask information after nodes." << std::endl
                      << "\tlandmask ne" << std::endl
		    << "where ne is the number of elements. Should match nx * ny" << std::endl
		    << "If all is land, just provide " << std::endl
		    << "   landmask 0" <<  std::endl;
            abort();
	}
	size_t ne;
	IN >> ne;

	if (ne == 0)
	  {
	    landmask.resize(nx*ny,true);
	  }
	else 
	  {
	    assert(ne == nx * ny);
	    landmask.resize(nx*ny,false);
	    bool lm;
	    for (size_t i=0;i<nx*ny;++i)
	      {
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
		      << "Then, for each edge we expect a line with two entries: id-of-element [0,1,2,3] for the edge" << std::endl;
            abort();
        }
        size_t nd;
        IN >> nd;

        if (statuslog > 0)
            std::cout << "reading " << nd << " dirichlet segments" << std::endl;

        for (size_t i = 0; i < nd; ++i) {
	  size_t n0, n1;
	  IN >> n0 >> n1; // read the element and the side
	  assert (n0<nx*ny);
	  assert (n0>=0);

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
                      << "where nd is the number of periodic edges" << std::endl
		      << "For each periodic edge: element-id1, element-id2, edge-of-element1, edge-of-element2" << std::endl;
            std::cerr << "Instead, got: " << status << std::endl;
            abort();
        }
        IN >> nd;
	if (statuslog > 0)
	  std::cout << "reading " << nd << " periodic edges" << std::endl;

        periodic.clear();
        periodic.resize(nd);
        for (size_t i = 0; i < nd; ++i) {
	    periodic[i].clear();
	    IN >> periodic[i].eid1 >> periodic[i].eid2 >> periodic[i].eoe1 >> periodic[i].eoe2;

	    size_t ix = periodic[i].eid1 % nx;
	    size_t iy = periodic[i].eid1 / nx;
	    
	    if      (periodic[i].eoe1 == 0) // bottom 
	      periodic[i].edgeid = nx*iy + ix;
	    else if (periodic[i].eoe1 == 2) // top
	      periodic[i].edgeid = nx*(iy+1) + ix;
	    else if (periodic[i].eoe1 == 3) // left 
	      periodic[i].edgeid = (nx+1)*iy + ix;
	    else if (periodic[i].eoe1 == 1) // right
	      periodic[i].edgeid = (nx+1)*iy + ix + 1;
	    
        }
    }

    IN.close();



    if (statuslog > 0) {
        std::cout << "ParametricMesh :: read mesh file " << fname << std::endl
                  << "             nx,ny = " << nx << " , " << ny << std::endl
                  << "             " << nelements << " elements,  " << nnodes << " nodes" << std::endl;
    }
}

/*!
 * Copy the coordinate arrays from the arguments.
 *
 * @param coord1 x in metres or longitude in radians
 * @param coord2 y in metres or latitude in radians
 */
void ParametricMesh::coordinatesFromModelArray(const ModelArray& coords)
{
    // Fill in the array sizes from the ModelArray dimensions
    nx = ModelArray::size(ModelArray::Dimension::X);
    ny = ModelArray::size(ModelArray::Dimension::Y);
    nelements = nx * ny;
    nnodes = (nx + 1) * (ny + 1);
    vertices.resize(nnodes, 2);
    for (size_t idx = 0; idx < nnodes; ++idx) {
        vertices(idx, 0) = coords.components(idx)[0];
        vertices(idx, 1) = coords.components(idx)[1];
    }
}

/*!
 * Copy the landmask from the passed ModelArray.
 *
 * @param mask the ModelArray containing the mask to be used.
 */
void ParametricMesh::landmaskFromModelArray(const ModelArray& mask)
{
    landmask.resize(nelements);
    for (size_t idx = 0; idx < mask.trueSize(); ++idx)
    {
        landmask[idx] = (mask[idx] == 1.);
    }
}

/*!
 * Add to the dirichlet arrays according to the stored landmask.
 */
void ParametricMesh::dirichletFromMask()
{
    // Edges are accessed in the order: BOTTOM, RIGHT, TOP, LEFT. See also ParametricMesh::edges.
    const std::array<size_t, N_EDGE> startX = {0, 0, 0, 1};
    const std::array<size_t, N_EDGE> stopX = {nx, nx - 1, nx, nx};
    const std::array<size_t, N_EDGE> startY = {1, 0, 0, 0};
    const std::array<size_t, N_EDGE> stopY = {ny, ny, ny - 1, ny};
    const std::array<int, N_EDGE> deltaIdx = {-static_cast<int>(nx), 1, static_cast<int>(nx), -1};

    // Loop over edges
    for (Edge edge : edges) {
        for (size_t j = startY[edge]; j < stopY[edge]; ++j) {
            for (size_t i = startX[edge]; i < stopX[edge]; ++i) {
                size_t idx = ModelArray::indexFromLocation(ModelArray::Type::H, {i, j});
                if (!landmask[idx]) continue;
                // mask(i, j) is ocean. Check the appropriate neighbour
                if (!landmask[idx + deltaIdx[edge]]) {
                    dirichlet[edge].push_back(idx);
                }
            }
        }
        sortDirichlet(edge);
    }
}

/*!
 * Add to the dirichlet arrays due to the domain edges according to an edge index.
 *
 * @param edge index of the edge to add closed boundary conditions to.
 */
void ParametricMesh::dirichletFromEdge(Edge edge)
{
    // BOTTOM, RIGHT, TOP, LEFT
    const std::array<size_t, N_EDGE> start = {0, nx - 1, nelements - nx, 0};
    const std::array<size_t, N_EDGE> stop = {nx, nelements, nelements, nelements};
    const std::array<size_t, N_EDGE> stride = {1, nx, 1, nx};

    for (size_t idx = start[edge]; idx < stop[edge]; idx += stride[edge]) {
        if (landmask[idx]) {
            dirichlet[edge].push_back(idx);
        }
    }
    sortDirichlet(edge);
}

/*!
 * Sort all the dirichlet arrays, so the element indices are ordered.
 */
void ParametricMesh::sortDirichlet()
{
    for (ParametricMesh::Edge edge : ParametricMesh::edges) {
        sortDirichlet(edge);
    }
}

/*!
 * Sort the dirichlet array of one particular edge.
 *
 * @param edge the edge to be sorted.
 */
void ParametricMesh::sortDirichlet(Edge edge)
{
    std::sort(dirichlet[edge].begin(), dirichlet[edge].end());
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
