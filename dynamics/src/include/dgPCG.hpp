/*!
 * @file dgPCG.hpp
 *
 * Matrix-free preconditioned conjugate gradient (PCG) method on DGVectors. The operator
 * A and the preconditioner P^-1 are passed in as functions. All vector operations are
 * local (cell by cell) and OpenMP-parallel.
 */

#ifndef __DGPCG_HPP
#define __DGPCG_HPP

#include "dgVector.hpp"

#include <cmath>

namespace Nextsim {

//! Global scalar product (x, y) = sum over cells c of x_c . y_c
template <int DG> FloatType dgDot(const DGVector<DG>& x, const DGVector<DG>& y)
{
    FloatType s = 0.0;
#pragma omp parallel for reduction(+ : s)
    for (size_t eid = 0; eid < (size_t)x.rows(); ++eid)
        s += x.row(eid).dot(y.row(eid));
    return s;
}

/*!
 * Solves A u = F by the preconditioned conjugate gradient method.
 *
 *   applyA(x, y)    : computes y = A x
 *   applyPrec(r, z) : computes z = P^-1 r
 *   u               : in: initial guess, out: solution
 *   tol             : stop when ||F - A u|| / ||F|| < tol
 *   maxit           : maximum number of iterations
 *   err             : out: final relative residual ||F - A u|| / ||F||
 *
 * Returns the number of iterations (maxit if not converged).
 */
template <int DG, class ApplyA, class ApplyPrec>
int PCG(const ApplyA& applyA, const ApplyPrec& applyPrec, const DGVector<DG>& F, DGVector<DG>& u,
    const FloatType tol, const int maxit, FloatType& err)
{
    const size_t ncells = F.rows();
    DGVector<DG> r, z, p, Ap;
    r.resize(ncells, DG);
    z.resize(ncells, DG);
    p.resize(ncells, DG);
    Ap.resize(ncells, DG);

    const FloatType normF = std::sqrt(dgDot(F, F));
    if (normF == 0.0) { // A u = 0  =>  u = 0
        u.setZero();
        err = 0.0;
        return 0;
    }

    // r = F - A u
    applyA(u, Ap);
#pragma omp parallel for
    for (size_t eid = 0; eid < ncells; ++eid)
        r.row(eid) = F.row(eid) - Ap.row(eid);

    err = std::sqrt(dgDot(r, r)) / normF;
    if (err < tol)
        return 0;

    // z = P^-1 r,  p = z
    applyPrec(r, z);
    p = z;
    FloatType rz = dgDot(r, z);

    for (int it = 1; it <= maxit; ++it) {
        applyA(p, Ap);
        const FloatType alpha = rz / dgDot(p, Ap);

        // u += alpha p,  r -= alpha A p
#pragma omp parallel for
        for (size_t eid = 0; eid < ncells; ++eid) {
            u.row(eid) += alpha * p.row(eid);
            r.row(eid) -= alpha * Ap.row(eid);
        }

        err = std::sqrt(dgDot(r, r)) / normF;
        if (err < tol)
            return it;

        applyPrec(r, z);
        const FloatType rz_new = dgDot(r, z);
        const FloatType beta = rz_new / rz;
        rz = rz_new;

        // p = z + beta p
#pragma omp parallel for
        for (size_t eid = 0; eid < ncells; ++eid)
            p.row(eid) = z.row(eid) + beta * p.row(eid);
    }
    return maxit;
}

} // namespace Nextsim

#endif /* __DGPCG_HPP */
