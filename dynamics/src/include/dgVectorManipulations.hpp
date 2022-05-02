/*!
 * @file dgVectorManipulations.hpp
 * @date 1 Mar 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#ifndef __DGVECTOR_MANIPULATIONS_HPP
#define __DGVECTOR_MANIPULATIONS_HPP

#include "dgVector.hpp"
#include <iomanip>

namespace Nextsim {
/**
* Average onto an edge between two adjacent elements 1/2 (c1 + c2)
* _X refers to a horizontal edge, _Y to a vertical edge
**/
LocalEdgeVector<1> average_to_X(const LocalCellVector<1>& down, const LocalCellVector<1>& up)
{
    return 0.5 * (down + up);
}

LocalEdgeVector<2> average_to_X(const LocalCellVector<3>& down, const LocalCellVector<3>& up)
{
    return LocalEdgeVector<2>(0.5 * (down(0, 0) + up(0, 0)) + 0.25 * (down(0, 2) - up(0, 2)),
        0.5 * (down(0, 1) + up(0, 1)));
}

LocalEdgeVector<3> average_to_X(const LocalCellVector<6>& down, const LocalCellVector<6>& up)
{
    return LocalEdgeVector<3>(0.5
            * (up(0, 0) - 0.5 * up(0, 2) + 1. / 6. * up(0, 4)
                + (down(0, 0) + 0.5 * down(0, 2) + 1. / 6. * down(0, 4))),
        0.5 * (up(0, 1) - 0.5 * up(0, 5) + (down(0, 1) + 0.5 * down(0, 5))),
        0.5 * (up(0, 3) + down(0, 3)));
}

LocalEdgeVector<1> average_to_Y(
    const Eigen::Matrix<double, 1, 1>& left, const Eigen::Matrix<double, 1, 1>& right)
{
    return 0.5 * (left + right);
}

LocalEdgeVector<2> average_to_Y(
    const Eigen::Matrix<double, 1, 3>& left, const Eigen::Matrix<double, 1, 3>& right)
{
    return LocalEdgeVector<2>(0.5 * (left(0, 0) + right(0, 0)) + 0.25 * (left(0, 1) - right(0, 1)),
        0.5 * (left(0, 2) + right(0, 2)));
}

LocalEdgeVector<3> average_to_Y(
    const Eigen::Matrix<double, 1, 6>& left, const Eigen::Matrix<double, 1, 6>& right)
{
    return LocalEdgeVector<3>(0.5
            * (right(0, 0) - 0.5 * right(0, 1) + 1. / 6. * right(0, 3)
                + (left(0, 0) + 0.5 * left(0, 1) + 1. / 6. * left(0, 3))),
        0.5 * (right(0, 1) - 0.5 * right(0, 5) + (left(0, 1) + 0.5 * left(0, 5))),
        0.5 * (right(0, 4) + left(0, 4)));
}

/**
 * Assemble the jump on an edge between two adjacent elements
 * [[u]] = u_1 n_1 + u_2 n_2
 * _X refers to a horizontal edge with n = (0, +/- 1)
 *    only the y-component is given, [u] = lower - upper
 * _Y refers to a vertical edge with n = (+/- 1,0)
 *    only the x-compoennt is given, [u] = left - right
 **/
template <int DG>
LocalEdgeVector<DG> jump_to_X(
    const LocalCellVector<DG>& down, const LocalCellVector<DG>& up);

template <>
LocalEdgeVector<1> jump_to_X(const LocalCellVector<1>& down, const LocalCellVector<1>& up)
{
    return down - up;
}
template <int DG>
LocalEdgeVector<DG> jump_to_Y(
    const LocalCellVector<DG>& left, const LocalCellVector<DG>& right);

template <>
LocalEdgeVector<1> jump_to_Y(const LocalCellVector<1>& left, const LocalCellVector<1>& right)
{
    return left - right;
}

//= -(phi(ic+1,0)-phi(ic,0));

//! Functions to average a DG cell vector onto the edges

void average_to_edges_Y(const Mesh& mesh, EdgeVector<1>& edgevector,
    const CellVector<1>& cellvector, bool settozero=true, bool settozeroonouteredges=false)
{
    assert(edgevector.edgetype == EdgeType::Y);

    if (settozero)
        edgevector.zero();
#pragma omp parallel for
    for (size_t iy = 0; iy < mesh.ny; ++iy) {
        size_t ey = iy * (mesh.nx + 1); // first edge index in row
        size_t cy = iy * mesh.nx; // first cell index in row
        edgevector.block(ey, 0, mesh.nx, 1) += 0.5 * cellvector.block(cy, 0, mesh.nx, 1);
        edgevector.block(ey + 1, 0, mesh.nx, 1) += 0.5 * cellvector.block(cy, 0, mesh.nx, 1);
    }

    if (settozeroonouteredges) { // outer edges are set to zero (Dirichlet-Values)
#pragma omp parallel for
        for (size_t iy = 0; iy < mesh.ny; ++iy) {
            edgevector(iy * (mesh.nx + 1), 0) = 0;
            edgevector(iy * (mesh.nx + 1) + mesh.nx, 0) = 0;
        }
    } else { // outer edges multiplied by 2
#pragma omp parallel for
        for (size_t iy = 0; iy < mesh.ny; ++iy) {
            edgevector(iy * (mesh.nx + 1), 0) *= 2.0;
            edgevector(iy * (mesh.nx + 1) + mesh.nx, 0) *= 2.0;
        }
    }
}

void average_to_edges_Y(const Mesh& mesh, EdgeVector<2>& edgevector,
    const CellVector<3>& cellvector, bool settozero=true, bool settozeroonouteredges=false)
{
    assert(edgevector.edgetype == EdgeType::Y);

    if (settozero)
        edgevector.zero();

#pragma omp parallel for
    for (size_t iy = 0; iy < mesh.ny; ++iy) {
        size_t ey = iy * (mesh.nx + 1); // first edge index in row
        size_t cy = iy * mesh.nx; // first cell index in row

        edgevector.block(ey, 0, mesh.nx, 1)
            += 0.5 * cellvector.block(cy, 0, mesh.nx, 1); // constant
        edgevector.block(ey + 1, 0, mesh.nx, 1) += 0.5 * cellvector.block(cy, 0, mesh.nx, 1);

        edgevector.block(ey, 0, mesh.nx, 1) += -0.25 * cellvector.block(cy, 1, mesh.nx, 1); // x-1/2
        edgevector.block(ey + 1, 0, mesh.nx, 1) += 0.25 * cellvector.block(cy, 1, mesh.nx, 1);

        edgevector.block(ey, 1, mesh.nx, 1) += 0.5 * cellvector.block(cy, 2, mesh.nx, 1); // x-1/2
        edgevector.block(ey + 1, 1, mesh.nx, 1) += 0.5 * cellvector.block(cy, 2, mesh.nx, 1);
    }

    if (settozeroonouteredges) {
#pragma omp parallel for
        for (size_t iy = 0; iy < mesh.ny; ++iy) {
            edgevector(iy * (mesh.nx + 1), 0) = 0;
            edgevector(iy * (mesh.nx + 1), 1) = 0;
            edgevector(iy * (mesh.nx + 1) + mesh.nx, 0) = 0;
            edgevector(iy * (mesh.nx + 1) + mesh.nx, 1) = 0;
        }
    } else {
#pragma omp parallel for
        for (size_t iy = 0; iy < mesh.ny; ++iy) {
            edgevector(iy * (mesh.nx + 1), 0) *= 2.0;
            edgevector(iy * (mesh.nx + 1), 1) *= 2.0;
            edgevector(iy * (mesh.nx + 1) + mesh.nx, 0) *= 2.0;
            edgevector(iy * (mesh.nx + 1) + mesh.nx, 1) *= 2.0;
        }
    }
}

void average_to_edges_Y(const Mesh& mesh, EdgeVector<3>& edgevector,
    const CellVector<6>& cellvector, bool settozero=true, bool settozeroonouteredges=false)
{
    assert(edgevector.edgetype == EdgeType::Y);

    if (settozero)
        edgevector.zero();
#pragma omp parallel for
    for (size_t iy = 0; iy < mesh.ny; ++iy) {
        size_t ic = iy * mesh.nx; // first index of left cell in row
        size_t ie = iy * (mesh.nx + 1); // first index of edge

        for (size_t i = 0; i < mesh.nx; ++i, ic += 1, ie += 1) // run over all elements
        {
            // add to left
            edgevector.block<1, 3>(ie, 0) += LocalEdgeVector<3>(
                0.5 * ((cellvector(ic, 0) - 0.5 * cellvector(ic, 1) + 1. / 6. * cellvector(ic, 3))),
                0.5 * ((cellvector(ic, 2) - 0.5 * cellvector(ic, 5))), 0.5 * (cellvector(ic, 4)));
            // add to right
            edgevector.block<1, 3>(ie + 1, 0) += LocalEdgeVector<3>(
                0.5 * ((cellvector(ic, 0) + 0.5 * cellvector(ic, 1) + 1. / 6. * cellvector(ic, 3))),
                0.5 * ((cellvector(ic, 2) + 0.5 * cellvector(ic, 5))), 0.5 * (cellvector(ic, 4)));
        }
    }

    if (settozeroonouteredges) {
#pragma omp parallel for
        for (size_t iy = 0; iy < mesh.ny; ++iy) {
            for (int c = 0; c < 3; ++c) {
                edgevector(iy * (mesh.nx + 1), c) = 0;
                edgevector(iy * (mesh.nx + 1) + mesh.nx, c) = 0;
            }
        }
    } else {
#pragma omp parallel for
        for (size_t iy = 0; iy < mesh.ny; ++iy) {
            for (int c = 0; c < 3; ++c) {
                edgevector(iy * (mesh.nx + 1), c) *= 2.0;
                edgevector(iy * (mesh.nx + 1) + mesh.nx, c) *= 2.0;
            }
        }
    }
}

void average_to_edges_X(const Mesh& mesh, EdgeVector<1>& edgevector,
    const CellVector<1>& cellvector, bool settozero=true, bool settozeroonouteredges=false)
{
    assert(edgevector.edgetype == EdgeType::X);

    if (settozero)
        edgevector.zero();
    // #pragma omp parallel for !!! overlap in writing. in each row we add to top
    // and bottom
    for (size_t iy = 0; iy < mesh.ny; ++iy) {
        size_t ex = iy * mesh.nx; // index of edge and cell x-vector
        edgevector.block(ex, 0, mesh.nx, 1) += 0.5 * cellvector.block(ex, 0, mesh.nx, 1);
        edgevector.block(ex + mesh.nx, 0, mesh.nx, 1) += 0.5 * cellvector.block(ex, 0, mesh.nx, 1);
    }

    if (settozeroonouteredges) {
#pragma omp parallel for
        for (size_t ix = 0; ix < mesh.nx; ++ix) {
            edgevector(ix, 0) = 0;
            edgevector(mesh.n + ix, 0) = 0;
        }
    } else {
#pragma omp parallel for
        for (size_t ix = 0; ix < mesh.nx; ++ix) {
            edgevector(ix, 0) *= 2.0;
            edgevector(mesh.n + ix, 0) *= 2.0;
        }
    }
}


void average_to_edges_X(const Mesh& mesh, EdgeVector<2>& edgevector,
    const CellVector<3>& cellvector, bool settozero=true, bool settozeroonouteredges=false)
{
    assert(edgevector.edgetype == EdgeType::X);

    if (settozero)
        edgevector.zero();
    //#pragma omp parallel for !!! see above!
    for (size_t iy = 0; iy < mesh.ny; ++iy) {
        size_t ex = iy * mesh.nx; // index of edge x-vector
        edgevector.block(ex, 0, mesh.nx, 1)
            += 0.5 * cellvector.block(ex, 0, mesh.nx, 1); // constants
        edgevector.block(ex + mesh.nx, 0, mesh.nx, 1) += 0.5 * cellvector.block(ex, 0, mesh.nx, 1);

        edgevector.block(ex, 1, mesh.nx, 1) += 0.5 * cellvector.block(ex, 1, mesh.nx, 1); // x-1/2
        edgevector.block(ex + mesh.nx, 1, mesh.nx, 1) += 0.5 * cellvector.block(ex, 1, mesh.nx, 1);

        edgevector.block(ex, 0, mesh.nx, 1) += -0.25 * cellvector.block(ex, 2, mesh.nx, 1); // y-1/2
        edgevector.block(ex + mesh.nx, 0, mesh.nx, 1) += 0.25 * cellvector.block(ex, 2, mesh.nx, 1);
    }

    if (settozeroonouteredges) {
#pragma omp parallel for
        for (size_t ix = 0; ix < mesh.nx; ++ix) {
            edgevector(ix, 0) = 0;
            edgevector(ix, 1) = 0;
            edgevector(mesh.n + ix, 0) = 0;
            edgevector(mesh.n + ix, 1) = 0;
        }
    } else {
#pragma omp parallel for
        for (size_t ix = 0; ix < mesh.nx; ++ix) {
            edgevector(ix, 0) *= 2.0;
            edgevector(ix, 1) *= 2.0;
            edgevector(mesh.n + ix, 0) *= 2.0;
            edgevector(mesh.n + ix, 1) *= 2.0;
        }
    }
}


void average_to_edges_X(const Mesh& mesh, EdgeVector<3>& edgevector,
    const CellVector<6>& cellvector, bool settozero=true, bool settozeroonouteredges=false)
{
    assert(edgevector.edgetype == EdgeType::X);

    if (settozero)
        edgevector.zero();

#pragma omp parallel for
    for (size_t ix = 0; ix < mesh.nx; ++ix) {
        size_t ic = ix; // first element
        size_t ie = ix; // first edge
        for (size_t i = 0; i < mesh.ny; ++i, ic += mesh.nx, ie += mesh.nx) {
            // add to bottom
            edgevector.block<1, 3>(ie, 0) += LocalEdgeVector<3>(
                0.5 * ((cellvector(ic, 0) - 0.5 * cellvector(ic, 2) + 1. / 6. * cellvector(ic, 4))),
                0.5 * ((cellvector(ic, 1) - 0.5 * cellvector(ic, 5))), 0.5 * (cellvector(ic, 3)));
            // add to top
            edgevector.block<1, 3>(ie + mesh.nx, 0) += LocalEdgeVector<3>(
                0.5 * ((cellvector(ic, 0) + 0.5 * cellvector(ic, 2) + 1. / 6. * cellvector(ic, 4))),
                0.5 * ((cellvector(ic, 1) + 0.5 * cellvector(ic, 5))), 0.5 * (cellvector(ic, 3)));
        }
    }
    if (settozeroonouteredges) {
#pragma omp parallel for
        for (size_t ix = 0; ix < mesh.nx; ++ix) {
            edgevector(ix, 0) = 0;
            edgevector(ix, 1) = 0;
            edgevector(ix, 2) = 0;
            edgevector(mesh.n + ix, 0) = 0;
            edgevector(mesh.n + ix, 1) = 0;
            edgevector(mesh.n + ix, 2) = 0;
        }
    } else {
#pragma omp parallel for
        for (size_t ix = 0; ix < mesh.nx; ++ix) {
            edgevector(ix, 0) *= 2.0;
            edgevector(ix, 1) *= 2.0;
            edgevector(ix, 2) *= 2.0;
            edgevector(mesh.n + ix, 0) *= 2.0;
            edgevector(mesh.n + ix, 1) *= 2.0;
            edgevector(mesh.n + ix, 2) *= 2.0;
        }
    }
}

//   template <int DG>
//   void jump_to_edges(const Mesh &mesh,
//                      EdgeVector<DG> &edgevector,
//                      const CellVector<DG> &cellvector,
//                      bool settozero = true,
//                      bool settozeroonouteredges = false);

//   template <>
//   void jump_to_edges(const Mesh &mesh,
//                      EdgeVector<0> &edgevector,
//                      const CellVector<0> &cellvector,
//                      bool settozero,
//                      bool settozeroonouteredges)
//   {
//     if (settozero)
//       edgevector.zero();
// #pragma omp parallel for
//     for (size_t iy = 0; iy < mesh.ny; ++iy)
//     {
//       size_t ex = iy * mesh.nx; // index of edge x-vector
//       edgevector.X.block(ex, 0, mesh.nx, 1) -=
//           cellvector.block(ex, 0, mesh.nx, 1);
//       edgevector.X.block(ex + mesh.nx, 0, mesh.nx, 1) +=
//           cellvector.block(ex, 0, mesh.nx, 1);
//       size_t ey = iy * (mesh.nx + 1);
//       edgevector.Y.block(ey, 0, mesh.nx, 1) -=
//           cellvector.block(ex, 0, mesh.nx, 1);
//       edgevector.Y.block(ey + 1, 0, mesh.nx, 1) +=
//           cellvector.block(ex, 0, mesh.nx, 1);
//     }

//     if (settozeroonouteredges)
//     {
// #pragma omp parallel for
//       for (size_t ix = 0; ix < mesh.nx; ++ix)
//       {
//         edgevector.X(ix, 0) = 0;
//         edgevector.X(mesh.n - mesh.nx + ix, 0) = 0;
//       }

// #pragma omp parallel for
//       for (size_t iy = 0; iy < mesh.ny; ++iy)
//       {
//         edgevector.Y(iy * (mesh.nx + 1), 0) = 0;
//         edgevector.Y(iy * (mesh.nx + 1) + mesh.nx, 0) = 0;
//       }
//     }
//   }
//   template <>
//   void jump_to_edges(const Mesh &mesh,
//                      EdgeVector<1> &edgevector,
//                      const CellVector<1> &cellvector,
//                      bool settozero,
//                      bool settozeroonouteredges)
//   {
//     if (settozero)
//       edgevector.zero();
// #pragma omp parallel for
//     for (size_t iy = 0; iy < mesh.ny; ++iy)
//     {
//       size_t ex = iy * mesh.nx; // index of edge x-vector
//       edgevector.X.block(ex, 0, mesh.nx, 1) +=
//           -cellvector.block(ex, 0, mesh.nx, 1); // constants
//       edgevector.X.block(ex + mesh.nx, 0, mesh.nx, 1) +=
//           cellvector.block(ex, 0, mesh.nx, 1);

//       edgevector.X.block(ex, 1, mesh.nx, 1) +=
//           -cellvector.block(ex, 1, mesh.nx, 1); // x-1/2
//       edgevector.X.block(ex + mesh.nx, 1, mesh.nx, 1) +=
//           cellvector.block(ex, 1, mesh.nx, 1);

//       edgevector.X.block(ex, 0, mesh.nx, 1) +=
//           0.5 * cellvector.block(ex, 2, mesh.nx, 1); // y-1/2
//       edgevector.X.block(ex + mesh.nx, 0, mesh.nx, 1) +=
//           -0.5 * cellvector.block(ex, 2, mesh.nx, 1);

//       size_t ey = iy * (mesh.nx + 1);
//       edgevector.Y.block(ey, 0, mesh.nx, 1) +=
//           -cellvector.block(ex, 0, mesh.nx, 1); // constant
//       edgevector.Y.block(ey + 1, 0, mesh.nx, 1) +=
//           cellvector.block(ex, 0, mesh.nx, 1);

//       edgevector.Y.block(ey, 0, mesh.nx, 1) +=
//           0.5 * cellvector.block(ex, 1, mesh.nx, 1); // x-1/2
//       edgevector.Y.block(ey + 1, 0, mesh.nx, 1) +=
//           -0.5 * cellvector.block(ex, 1, mesh.nx, 1);

//       edgevector.Y.block(ey, 1, mesh.nx, 1) +=
//           -cellvector.block(ex, 2, mesh.nx, 1); // x-1/2
//       edgevector.Y.block(ey + 1, 1, mesh.nx, 1) +=
//           cellvector.block(ex, 2, mesh.nx, 1);
//     }

//     if (settozeroonouteredges)
//     {
// #pragma omp parallel for
//       for (size_t ix = 0; ix < mesh.nx; ++ix)
//       {
//         edgevector.X(ix, 0) = 0;
//         edgevector.X(ix, 1) = 0;
//         edgevector.X(mesh.n + ix, 0) = 0;
//         edgevector.X(mesh.n + ix, 1) = 0;
//       }

// #pragma omp parallel for
//       for (size_t iy = 0; iy < mesh.ny; ++iy)
//       {
//         edgevector.Y(iy * (mesh.nx + 1), 0) = 0;
//         edgevector.Y(iy * (mesh.nx + 1), 1) = 0;
//         edgevector.Y(iy * (mesh.nx + 1) + mesh.nx, 0) = 0;
//         edgevector.Y(iy * (mesh.nx + 1) + mesh.nx, 1) = 0;
//       }
//     }
//   }

//   template <>
//   void jump_to_edges(const Mesh &mesh,
//                      EdgeVector<2> &edgevector,
//                      const CellVector<2> &cellvector,
//                      bool settozero,
//                      bool settozeroonouteredges)
//   {
//     abort();
//   }

//   template <int DG>
//   void print_cell_vector(const Mesh &mesh, const CellVector<DG>
//   &cellvector); template <> void print_cell_vector(const Mesh &mesh, const
//   CellVector<0> &cellvector)
//   {
//     for (int y = 0; y < mesh.ny; ++y)
//     {
//       for (int x = 0; x < mesh.nx; ++x)
//         std::cout << std::setw(6)
//                   << cellvector((mesh.ny - y - 1) * mesh.nx + x, 0);
//       std::cout << std::endl;
//     }
//   }

//   template <int DG>
//   void print_edge_vector(const Mesh &mesh, const EdgeVector<DG>
//   &edgevector); template <> void print_edge_vector(const Mesh &mesh, const
//   EdgeVector<0> &edgevector)
//   {
//     std::cout << "X:" << std::endl;
//     for (int y = 0; y <= mesh.ny; ++y)
//     {
//       for (int x = 0; x < mesh.nx; ++x)
//         std::cout << std::setw(6)
//                   << edgevector.X((mesh.ny - y) * mesh.nx + x, 0);
//       std::cout << std::endl;
//     }

//     std::cout << "Y:" << std::endl;
//     for (int y = 0; y < mesh.ny; ++y)
//     {
//       for (int x = 0; x <= mesh.nx; ++x)
//         std::cout << std::setw(6)
//                   << edgevector.X((mesh.ny - y - 1) * mesh.nx + x, 0);
//       std::cout << std::endl;
//     }
//   }

} /* namespace Nextsim */

#endif /* __DGVECTOR_MANIPULATIONS_HPP */
