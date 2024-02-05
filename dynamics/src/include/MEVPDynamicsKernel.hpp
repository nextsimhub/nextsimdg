/*!
 * @file MEVPDynamicsKernel.hpp
 *
 * @date 17 Feb 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Piotr Minakowski <piotr.minakowski@ovgu.de>
 */

#ifndef MEVPDYNAMICSKERNEL_HPP
#define MEVPDYNAMICSKERNEL_HPP

#include "VPCGDynamicsKernel.hpp"

namespace Nextsim {

template <int DGadvection> class MEVPDynamicsKernel : public VPCGDynamicsKernel<DGadvection> {
public:
    MEVPDynamicsKernel(ParametricMomentumMap<CGdegree>* pmapIn, const DynamicsParameters& paramsIn)
        : VPCGDynamicsKernel<DGadvection>(pmapIn, MEVPStressStep, paramsIn)
    {
    }

private:
    //! Rheology-Parameters
    Nextsim::VPParameters VP;
    MEVPStressUpdateStep MEVPStressStep;
    //! MEVP parameters
    double alpha = 1500.0;
    double beta = 1500.0;
};

} /* namespace Nextsim */

#endif /* MEVPDYNAMICSKERNEL_HPP */
