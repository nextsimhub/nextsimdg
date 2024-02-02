/*!
 * @file StressUpdateStep.hpp
 *
 * @date Feb 1, 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef STRESSUPDATESTEP_HPP
#define STRESSUPDATESTEP_HPP

#include <map>
#include <string>

#include "dgVector.hpp"
#include "DynamicsParameters.hpp"
#include "ParametricMesh.hpp"

namespace Nextsim {

template <int DGadvection, int DGstress>
class StressUpdateStep {
public:
    typedef std::array<DGVector<DGstress>&, N_TENSOR_ELEMENTS> SymmetricTensorVector;

    static const int nGauss = ( ((DGstress == 8) || (DGstress == 6) ) ? 3 : (DGstress == 3 ? 2 : -1));
    StressUpdateStep() = default;
    virtual ~StressUpdateStep() = default;
    virtual void stressUpdateHighOrder(const DynamicsParameters& params,
            const ParametricMesh& smesh,
            const SymmetricTensorVector& stress, const SymmetricTensorVector& strain,
            const DGVector<DGadvection>& h, const DGVector<DGadvection>& a,
            const double deltaT) = 0;
};

} /* namespace Nextsim */

#endif /* STRESSUPDATESTEP_HPP */
