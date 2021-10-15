/*----------------------------   mesh.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __mesh_H
#define __mesh_H
/*----------------------------   mesh.h     ---------------------------*/

#include <array>

namespace Nextsim {

typedef std::array<double, 2> Vertex;

class Mesh {
public:
    size_t nx, ny; // no of elements in x- and y-direction
    size_t n; // total number of elements

    double h; // spatial mesh size

    Mesh()
        : nx(0)
        , ny(0)
        , n(0)
        , h(0)
    {
    }

    Mesh(size_t NX, size_t NY, double H)
        : nx(NX)
        , ny(NY)
        , n(NX * NY)
        , h(H)
    {
        assert(nx > 0);
        assert(ny > 0);
        assert(n == nx * ny);
        assert(h > 0);
    }

    Vertex vertex(size_t ix, size_t iy) const
    {
        return { h * (ix + 0.5), h * (iy + 0.5) };
    }
};

} // namespace Nextsim

/*----------------------------   mesh.h     ---------------------------*/
/* end of #ifndef __mesh_H */
#endif
/*----------------------------   mesh.h     ---------------------------*/
