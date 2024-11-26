/*!
 * @file FreeDriftDynamicsKernel.hpp
 *
 * Implementation of "classic free drift", where we ignore all \rho h terms in the momentum
 * equation. This is equivalent to assuming that the ice is very thin.
 *
 * @date 19 Nov 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#ifndef FREEDRIFTDYNAMICSKERNEL_HPP
#define FREEDRIFTDYNAMICSKERNEL_HPP

#include "CGDynamicsKernel.hpp"

namespace Nextsim {

template <int DGadvection> class FreeDriftDynamicsKernel : public CGDynamicsKernel<DGadvection> {
    using DynamicsKernel<DGadvection, DGstressComp>::nSteps;
    using DynamicsKernel<DGadvection, DGstressComp>::hice;
    using DynamicsKernel<DGadvection, DGstressComp>::cice;
    using DynamicsKernel<DGadvection, DGstressComp>::advectionAndLimits;
    using DynamicsKernel<DGadvection, DGstressComp>::dgtransport;

    using CGDynamicsKernel<DGadvection>::u;
    using CGDynamicsKernel<DGadvection>::v;
    using CGDynamicsKernel<DGadvection>::uOcean;
    using CGDynamicsKernel<DGadvection>::vOcean;
    using CGDynamicsKernel<DGadvection>::uAtmos;
    using CGDynamicsKernel<DGadvection>::vAtmos;
    using CGDynamicsKernel<DGadvection>::applyBoundaries;

public:
    FreeDriftDynamicsKernel(const DynamicsParameters& paramsIn)
        : CGDynamicsKernel<DGadvection>()
        , params(paramsIn)
    {
    }

    virtual ~FreeDriftDynamicsKernel() = default;

    void update(const TimestepTime& tst) override
    {
        updateMomentum(tst);
        applyBoundaries();

        // Let DynamicsKernel handle the advection step
        advectionAndLimits(tst);
    };

protected:
    const DynamicsParameters& params;

    const double cosOceanAngle = cos(radians * params.oceanTurningAngle);
    const double sinOceanAngle = sin(radians * params.oceanTurningAngle);
    const double FOcean = params.COcean * params.rhoOcean;
    const double FAtm = params.CAtm * params.rhoAtm;
    const double NansenNumber = sqrt(FAtm / FOcean);

    void updateMomentum(const TimestepTime& tst) override
    {
#pragma omp parallel for
        for (int i = 0; i < u.rows(); ++i) {
            // Free drift ice velocity
            u(i) = uOcean(i)
                + NansenNumber * (uAtmos(i) * cosOceanAngle - vAtmos(i) * sinOceanAngle);
            v(i) = vOcean(i)
                + NansenNumber * (-uAtmos(i) * sinOceanAngle + vAtmos(i) * cosOceanAngle);
        }
    }
};
} /* namespace Nextsim */

#endif /* FREEDRIFTDYNAMICSKERNEL_HPP */
