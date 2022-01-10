/*----------------------------   cg2dg.hpp     ---------------------------*/
#ifndef __cg2dg_HPP
#define __cg2dg_HPP
/*----------------------------   cg2dg.hpp     ---------------------------*/

/*!
 * @file    cg2dg.hpp
 * @author  Thomas Richter <thomas.richter@ovgu.de>
 * @brief   Projecting a CG vector to a DG
 */

#include "cgvector.hpp"
#include "codegeneration_cg_to_dg.hpp"
#include "dgvector.hpp"

namespace Nextsim {

template <int CG, int DG>
void ProjectCG2DG(const Mesh& mesh, CellVector<DG>& dg, const CGVector<CG>& cg);

template <>
void ProjectCG2DG(const Mesh& mesh, CellVector<1>& dg, const CGVector<2>& cg)
{
    assert(static_cast<long int>((2 * mesh.nx + 1) * (2 * mesh.ny + 1)) == cg.rows());
    assert(static_cast<long int>(mesh.nx * mesh.ny) == dg.rows());

    const int cgshift = 2 * mesh.nx + 1; //!< Index shift for each row

    // parallelize over the rows
#pragma omp parallel for
    for (size_t row = 0; row < mesh.ny; ++row) {
        int dgi = mesh.nx * row; //!< Index of dg vector
        int cgi = 2 * cgshift * row; //!< Lower left index of cg vector

        for (size_t col = 0; col < mesh.nx; ++col, ++dgi, cgi += 2) {
            Eigen::Vector<double, 9> cg_local = {
                cg(cgi), cg(cgi + 1), cg(cgi + 2),
                cg(cgi + cgshift), cg(cgi + 1 + cgshift), cg(cgi + 2 + cgshift),
                cg(cgi + 2 * cgshift), cg(cgi + 1 + 2 * cgshift), cg(cgi + 2 + 2 * cgshift)
            };
            dg.row(dgi) = CG2_to_DG1 * cg_local;
        }
    }
}
template <>
void ProjectCG2DG(const Mesh& mesh, CellVector<2>& dg, const CGVector<2>& cg)
{
    assert(static_cast<long int>((2 * mesh.nx + 1) * (2 * mesh.ny + 1)) == cg.rows());
    assert(static_cast<long int>(mesh.nx * mesh.ny) == dg.rows());

    const int cgshift = 2 * mesh.nx + 1; //!< Index shift for each row

    // parallelize over the rows
#pragma omp parallel for
    for (size_t row = 0; row < mesh.ny; ++row) {
        int dgi = mesh.nx * row; //!< Index of dg vector
        int cgi = 2 * cgshift * row; //!< Lower left index of cg vector

        for (size_t col = 0; col < mesh.nx; ++col, ++dgi, cgi += 2) {
            Eigen::Vector<double, 9> cg_local = {
                cg(cgi), cg(cgi + 1), cg(cgi + 2),
                cg(cgi + cgshift), cg(cgi + 1 + cgshift), cg(cgi + 2 + cgshift),
                cg(cgi + 2 * cgshift), cg(cgi + 1 + 2 * cgshift), cg(cgi + 2 + 2 * cgshift)
            };

            dg.row(dgi) = CG2_to_DG2 * cg_local;
        }
    }
}

}

/*----------------------------   cg2dg.hpp     ---------------------------*/
/* end of #ifndef __cg2dg_HPP */
#endif
/*----------------------------   cg2dg.hpp     ---------------------------*/
