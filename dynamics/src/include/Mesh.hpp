/*!
 * @file Mesh.hpp
 * @date 1 Mar 2022
 * @author Thomas Richter <thomas.richter@ovgu.no>
 */

#ifndef __MESH_HPP
#define __MESH_HPP

#include <array>
#include <cassert>

namespace Nextsim {

typedef std::array<double, 2> Vertex;

class Mesh {
public:
    size_t nx, ny; // no of elements in x- and y-direction
    size_t n; // total number of elements

    double hx, hy; // spatial mesh size

    Mesh()
        : nx(0)
        , ny(0)
        , n(0)
        , hx(0)
        , hy(0)
    {
    }

    Mesh(size_t NX, size_t NY, double HX, double HY)
        : nx(NX)
        , ny(NY)
        , n(NX * NY)
        , hx(HX)
        , hy(HY)
    {
        assert(nx > 0);
        assert(ny > 0);
        assert(n == nx * ny);
        assert(hx > 0);
        assert(hy > 0);
    }

    void BasicInit(size_t NX, size_t NY, double HX, double HY)
    {
        nx = NX;
        ny = NY;
        hx = HX;
        hy = HY;
        n = NX * NY;
    }

    Vertex midpoint(size_t ix, size_t iy) const { return { hx * (ix + 0.5), hy * (iy + 0.5) }; }
};

} /* namespace Nextsim */

#endif /* __MESH_HPP */
