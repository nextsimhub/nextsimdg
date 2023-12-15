/*!
 * @file ParametricMesh.hpp
 * @date 9 Juli 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#ifndef __PARAMETRICMESH_HPP
#define __PARAMETRICMESH_HPP

#include <Eigen/Dense>
#include <array>
#include <cassert>
#include <cstddef>

#include "include/ModelArray.hpp"
#include "NextsimDynamics.hpp"
#include <iostream>
#include <string>
#include <vector>

namespace Nextsim {

inline constexpr double SQR(double x) { return x * x; }

/*!
 * Stores the spatial mesh of the domain.
 *
 * The mesh has a structured topology with nelements = nx * ny elements
 * and nnodes = (nx+1) * (ny+1) nodes.
 * Everything is sorted from lower left to upper right, using
 *   cellindex = ny     * rowindex + colindex
 *   nodeindex = (ny+1) * rowindex + colindex
 *
 * Coordinates of the nnodes vertices are stored as
 *
 *   Eigen::Matrix< NextSim::FloatType , nnodes, 2 > vertices;
 *
 */

typedef std::array<double, 2> Vertex;

class ParametricMesh {
public:
    int statuslog; //!< -1 no output, 1 full status output

    enum Edge {
        BOTTOM,
        RIGHT,
        TOP,
        LEFT,
        N_EDGE
    };

    static constexpr std::array<Edge, N_EDGE> edges = {BOTTOM, RIGHT, TOP, LEFT};

  COORDINATES CoordinateSystem; //! CARTESIAN or SPHERICAL
  
    size_t nx, ny; //!< no of elements in x- and y-direction
    size_t nnodes; //!< total number of nodes
    size_t nelements; //!< total number of nodes

    Eigen::Matrix<Nextsim::FloatType, Eigen::Dynamic, 2> vertices; // stores the

    /*!
     * Stores the Dirichlet boundary information
     *
     * Dirichlet-Information is stored in 4 vectors
     * dirichlet[0]: List of elements that have a Dirichlet-Boundary on the bottom
     * dirichlet[1]: List of elements that have a Dirichlet-Boundary on the right
     * dirichlet[2]: List of elements that have a Dirichlet-Boundary on the top
     * dirichlet[3]: List of elements that have a Dirichlet-Boundary on the left
     *
     * Each of the vectors is processed in parallel
     */
  std::array<std::vector<size_t>, 4> dirichlet;
  
  /*! 
     * Periodic:
     *   [0]=0 for X-edge (bottom/top) and [0]=1 for Y-edge (left/right)
     *   [1] the element on the left / bottom
     *   [2] the element on the right / top
     *   [3] the edge id
     */
  std::vector<std::vector<std::array<size_t, 4>>> periodic;

  /*!
   * Landmask for the elements. true for land, false ice
   * if no landmask is given, all  elements are ice
   */
  std::vector<bool> landmask;

  ParametricMesh(const COORDINATES coords, int loglevel = -1)
    : CoordinateSystem (coords)
    , statuslog(loglevel)
        , nx(-1)
        , ny(-1)
        , nnodes(-1)
        , nelements(-1)
    {
    }

  /*! 
   * resets mesh to initial state
   */
  void reset()
  {
    nx = -1;
    ny = -1;
    nnodes = -1;
    nelements = -1;
    for (size_t i=0;i<4;++i)
      dirichlet[i].clear();
    periodic.clear();
    landmask.clear();
  }
  
  /*!
    * Reads mesh from a .smesh - file
    *
    * File format:
    *
    * ParametricMesh 1.0 % Identifier & Version
    *
    * nx ny          % number of elements in x- and y- direction  << Version: 1.0
    * x1       y1    % x- and y- coordinates of vertices. first vertex lower left
    * x2       y2    % then to the right, then 2nd row, ... up to upper right
    * ...
    * xxnnodes ynnodes
    */
  void readmesh(std::string fname);


  /*!
   * changes from [-180,180] to [-pi,pi] and [-90,90] to [-pi/2,pi/2]
   * NO CHECKING IF REQUIRED!
   */
  void TransformToRadians()
  {
    vertices *= M_PI/180.0;      
  }

  /*!
   * Rotates the mesh such that the singularities are in Greenland / Antarctica at
   * 75N / 40W and 75S / 140E
   */
  void RotatePoleToGreenland()
  {
    for (size_t i=0;i<nnodes;++i)
      {
	const double x = cos(vertices(i,1))*cos(vertices(i,0));
	const double y = cos(vertices(i,1))*sin(vertices(i,0));
	const double z = sin(vertices(i,1));

	double aw =  40.0*M_PI/180.0;
	const double x1 = cos(aw)*x - sin(aw)*y;
	const double y1 = sin(aw)*x + cos(aw)*y;
	const double z1 = z;

	double bw =  15.0*M_PI/180.0;
	const double x2 = cos(bw)*x1 - sin(bw)*z1;
	const double y2 = y1;
	const double z2 = sin(bw)*x1 + cos(bw)*z1;
	
	vertices(i,1) = asin(z2);
	vertices(i,0) = atan2(y2,x2);
      }
  }
  //! Rotation back to normal
  void RotatePoleFromGreenland()
  {
    for (size_t i=0;i<nnodes;++i)
      {
	const double x = cos(vertices(i,1))*cos(vertices(i,0));
	const double y = cos(vertices(i,1))*sin(vertices(i,0));
	const double z = sin(vertices(i,1));
	
	double aw = -40.0*M_PI/180.0;
	double bw = -15*M_PI/180.0;

	const double x1 = cos(bw)*x - sin(bw)*z;
	const double y1 = y;
	const double z1 = sin(bw)*x + cos(bw)*z;

	const double x2 = cos(aw)*x1 - sin(aw)*y1;
	const double y2 = sin(aw)*x1 + cos(aw)*y1;
	const double z2 = z1;

	
	vertices(i,1) = asin(z2);
	vertices(i,0) = atan2(y2,x2);
      }
  }


  
    /*!
     * Returns index of lower left node in element with index eid
     */
    size_t eid2nid(const size_t eid) const
    {
        return static_cast<size_t>(eid / nx) * (nx + 1) + (eid % nx);
    }


  /*!
   * makes sure that the longitudes in the coordinate vector coords
   * do not jump form Pi to -Pi degrees. If a jump is detected, 
   * all cordinates are shifted to the range 0,3 Pi by correcting
   * negative ones
   */
  template<int N>
  void correctlongitude(Eigen::Matrix<Nextsim::FloatType, N, 2>& coords) const
  {
    assert(CoordinateSystem == SPHERICAL);
    bool problem = false;
    for (size_t i=1;i<N;++i)
      if (fabs(coords(0,0)-coords(i,0))>2.0/3.0 * M_PI)
	{
	  problem = true;
	  break;
	}
    if (problem)
      {
	for (size_t i=0;i<N;++i)
	  if (coords(i,0)<0)
	    coords(i,0) += 2.0 * M_PI;
      }
  }

  /*!
     * returns the coordinates of one element with index eid as 4 x 2 matrix
     *
     * If correctlongitude is set to true, the congitude coordinates will
     * be such that there is no jump (from Pi to -Pi) within the element
     */
  const Eigen::Matrix<Nextsim::FloatType, 4, 2> coordinatesOfElement(const size_t eid) const
    {
        const size_t nid = eid2nid(eid);
	assert(nid<vertices.rows());
	assert(nid+nx+2<vertices.rows());
	Eigen::Matrix<Nextsim::FloatType, 4, 2> coords = 
	  Eigen::Matrix<Nextsim::FloatType, 4, 2>
	  ({ { vertices(nid, 0),          vertices(nid, 1) },
	     { vertices(nid + 1, 0),      vertices(nid + 1, 1) },
	     { vertices(nid + nx + 1, 0), vertices(nid + nx + 1, 1) },
	     { vertices(nid + nx + 2, 0), vertices(nid + nx + 2, 1) } }
	    );
	if (CoordinateSystem == SPHERICAL)
	  correctlongitude(coords);

	return coords;
    }


    /*!
     * returns the coordinates of one edge with index eid as 2 x 2 matrix 
     * (first index is the index of the coordinate, second the x/y - value)
     *
     * If correctlongitude is set to true, the congitude coordinates will
     * be such that there is no jump (from Pi/2 to -Pi/2) within the element
     */
  const Eigen::Matrix<Nextsim::FloatType, 2, 2> coordinatesOfEdgeX(const size_t eid) const
    {
      const size_t ex = eid%nx;        //! x-index of the corresponding element
      const size_t ey = eid/nx;        //! y-index of the corresponding element 
      const size_t nid = ey*(nx+1)+ex; //! index of the node 

      Eigen::Matrix<Nextsim::FloatType, 2, 2> ecoords =
	Eigen::Matrix<Nextsim::FloatType, 2, 2> 
	({ { vertices(nid, 0),          vertices(nid, 1) },
	   { vertices(nid + 1, 0),      vertices(nid + 1, 1) } });

	if (CoordinateSystem == SPHERICAL)
	correctlongitude(ecoords);
      
      return ecoords;
    }
  const Eigen::Matrix<Nextsim::FloatType, 2, 2> coordinatesOfEdgeY(const size_t eid) const
    {
      const size_t ex = eid%(nx+1); //! x-index of the corresponding element (possibly nx+1)
      const size_t ey = eid/(nx+1); //! y-index of the corresponding element
      assert(ey<ny);
      const size_t nid = ey*(nx+1) + ex;

      Eigen::Matrix<Nextsim::FloatType, 2, 2> ecoords =
	Eigen::Matrix<Nextsim::FloatType, 2, 2> 
	({ { vertices(nid, 0),         vertices(nid, 1) },
	   { vertices(nid + nx+1, 0),  vertices(nid + nx+1, 1) } });

	if (CoordinateSystem == SPHERICAL)
	correctlongitude(ecoords);
      
      return ecoords;

    }
  
    /*!
     * return the area of the mesh element with index eid
     */
    double area(const size_t eid) const
    {
        const size_t nid = eid2nid(eid);

        const double a = (vertices.block<1, 2>(nid, 0) - vertices.block<1, 2>(nid + 1, 0)).squaredNorm(); // lower
        const double b = (vertices.block<1, 2>(nid + 1, 0) - vertices.block<1, 2>(nid + nx + 2, 0)).squaredNorm(); // right
        const double c = (vertices.block<1, 2>(nid + 1 + nx, 0) - vertices.block<1, 2>(nid + 2 + nx, 0)).squaredNorm(); // top
        const double d = (vertices.block<1, 2>(nid, 0) - vertices.block<1, 2>(nid + nx + 1, 0)).squaredNorm(); // left
        const double e = (vertices.block<1, 2>(nid, 0) - vertices.block<1, 2>(nid + nx + 2, 0)).squaredNorm(); // diag 1
        const double f = (vertices.block<1, 2>(nid + 1, 0) - vertices.block<1, 2>(nid + nx + 1, 0)).squaredNorm(); // diag 2

        return 0.25 * sqrt(4.0 * e * f - SQR(b + d - a - c));
    }
    /*!
     * return the mesh size (as square root of area)  of the mesh element eid
     */
    double h(const size_t eid) const
    {
        return sqrt(area(eid));
    }

    /*!
     * returns the cooridnate of a vertex by an index (ix,iy)
     * For CG1 this is directly the corresponding entry ( (nx+1) iy + ix ) in the vertices-array
     * For CG2 we linearly interpolate.
     */
    template <int CGdegree>
    Vertex coordinate(size_t ix, size_t iy) const
    {
        if (CGdegree == 1) {
            assert(ix < nx + 1);
            assert(iy < ny + 1);
            const size_t ii = (nx + 1) * iy + ix;
            return Vertex({ vertices(ii, 0), vertices(ii, 1) });
        } else if (CGdegree == 2) {
            assert(ix < 2 * nx + 1);
            assert(iy < 2 * ny + 1);

	    const size_t cx = ix / 2; // index of the element where the cooridnate belongs to
            const size_t cy = iy / 2;

	    const size_t mx = ix % 2; // where are we in the element?
	    const size_t my = iy % 2;

	    if ( (cx<nx) && (cy<ny) )
	      {
		const Eigen::Matrix<Nextsim::FloatType, 4, 2> coe = coordinatesOfElement(nx*cy+cx);
		
		if ((my == 0) && (mx == 0))
		  return Vertex({ coe(0,0),coe(0,1) });
		else if ((my == 0) && (mx == 1))
		  return Vertex({ 0.5*(coe(0,0)+coe(1,0)), 0.5*(coe(0,1)+coe(1,1)) });
		else if ((my == 1) && (mx == 0))
		  return Vertex({ 0.5*(coe(0,0)+coe(2,0)), 0.5*(coe(0,1)+coe(2,1)) });
		else if ((my == 1) && (mx == 1))
		  return Vertex({ 0.25*(coe(0,0)+coe(1,0)+coe(2,0)+coe(3,0)), 0.25*(coe(0,1)+coe(1,1)+coe(2,1)+coe(3,1)) });
		else
		  abort();
	      }
	    else if (cx<nx)
	      // the node is on the upper boundary of  the mesh
	      // we take the element of the last row and access the upper nodes there
	      // it must be my = 0.
	      {
		assert(my == 0);
		const Eigen::Matrix<Nextsim::FloatType, 4, 2> coe = coordinatesOfElement(nx*(cy-1)+cx);
		
		if (mx == 0)
		  return Vertex({ coe(2,0),coe(2,1) });
		else if (mx == 1)
		  return Vertex({ 0.5*(coe(2,0)+coe(3,0)), 0.5*(coe(2,1)+coe(3,1)) });
		abort();
	      }
	    else if (cy<ny)
	      // the node is on the right boundary of  the mesh
	      // we take the element of the last col and access the right nodes there
	      // it must be mx = 0.
	      {
		assert(mx == 0);

		const Eigen::Matrix<Nextsim::FloatType, 4, 2> coe = coordinatesOfElement(nx*cy+cx-1);

		if (my == 0)
		  return Vertex({ coe(1,0),coe(1,1) });
		else if (my == 1)
		  return Vertex({ 0.5*(coe(1,0)+coe(3,0)),0.5*(coe(1,1)+coe(3,1)) });
		abort();
	      }
	    else // top right vertex
	      {
		assert( (cx == nx) && (ny == ny) );
		return Vertex( {vertices(nnodes-1,0), vertices(nnodes-1,1)} );
	      }
	    abort();
  
        } else
            abort();
    }

    /*!
     * computes the vector between two points n1, n2 (not normed)
     */
    Eigen::Matrix<Nextsim::FloatType, 1, 2> edgevector(const size_t n1, const size_t n2) const
  {
      //! In spherical coordinates (and greenland) we must check for the Pi -> -Pi jump
      if (CoordinateSystem == SPHERICAL)
	{
	  Eigen::Matrix<Nextsim::FloatType, 1, 2> dv = vertices.block<1,2>(n2,0) - vertices.block<1,2>(n1,0);

	  if (dv(0,0)>0.5*M_PI)
	    dv(0,0) -= 2.0 * M_PI;
	  if (dv(0,0)<-0.5*M_PI)
	    dv(0,0) += 2.0 * M_PI;
	  return dv;
	}
      else if (CoordinateSystem == CARTESIAN)
	return vertices.block<1,2>(n2,0) - vertices.block<1,2>(n1,0);
      else abort();      
    }

    /*!
     * Copy the coordinate arrays from the arguments.
     *
     * @param coord1 x in metres or longitude in radians
     * @param coord2 y in metres or latitude in radians
     */
    void coordinatesFromModelArray(const ModelArray& coord1, const ModelArray& coord2);

    /*!
     * Copy the landmask from the passed ModelArray
     */
    void landmaskFromModelArray(const ModelArray& mask);

    /*!
     * Add to the dirichlet arrays given a landmask.
     */
    void dirichletFromMask(const ModelArray& mask);

    /*!
     * Add to the dirichlet arrays according to an edge index.
     */
    void dirichletFromEdge(const ModelArray& mask, Edge edge);

    /*!
     * Sort the dirichlet arrays, so the indices are ordered.
     */
    void sortDirichlet();

    // Global access functions

    double hmin() const; //! returns the minimum mesh size
    double area() const; //! returns the area of the domain
};

} /* namespace Nextsim */

#endif /* __MESH_HPP */
