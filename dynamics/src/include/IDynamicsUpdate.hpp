/*!
 * @file IDynamicsUpdate.hpp
 *
 * @date Jan 19, 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IDYNAMICSUPDATE_HPP
#define IDYNAMICSUPDATE_HPP

#include "DynamicsParameters.hpp"
#include "ParametricMesh.hpp"
#include "dgVector.hpp"

#include <array>

namespace Nextsim {

class IDynamicsUpdate {
protected:
    //    static const size_t nSymMatrixEl = 3;
public:
    enum SymMatrixEl { _11, _12, _22, nSymMatrixEl };
    template <int CG, int DGstress, int DGadvection>
    void stressUpdateHighOrder(std::array<DGVector<DGstress>, nSymMatrixEl> S,
        const std::array<DGVector<DGstress>, nSymMatrixEl> E, const DGVector<DGadvection>& A,
        DGVector<DGadvection>& D, const double dt_mom)
        = 0;

protected:
    inline constexpr double SQR(double x) { return x * x; }
    const DynamicsParameters& params;
    const ParametricMesh& smesh;
};

}

#endif /* IDYNAMICSUPDATE_HPP */
