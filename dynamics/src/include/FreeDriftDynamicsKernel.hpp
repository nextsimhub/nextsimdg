/*!
 * @file FreeDriftDynamicsKernel.hpp
 *
 * Implementation of "classic free drift", where we ignore all \rho h terms in the momentum
 * equation. This is equivalent to assuming that the ice is very thin.
 *
 * @date 17 Feb 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#ifndef FREEDRIFTDYNAMICSKERNEL_HPP
#define FREEDRIFTDYNAMICSKERNEL_HPP

#include "CGDynamicsKernel.hpp"

namespace Nextsim {

template <int DGadvection> class FreeDriftDynamicsKernel : public CGDynamicsKernel<DGadvection> {
    using DynamicsKernel<DGadvection, DGstressDegree>::nSteps;
    using DynamicsKernel<DGadvection, DGstressDegree>::hice;
    using DynamicsKernel<DGadvection, DGstressDegree>::cice;
    using DynamicsKernel<DGadvection, DGstressDegree>::advectionAndLimits;
    using DynamicsKernel<DGadvection, DGstressDegree>::dgtransport;

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

    const double cosOceanAngle = cos(radians * params.ocean_turning_angle);
    const double sinOceanAngle = sin(radians * params.ocean_turning_angle);
    const double NansenNumber = sqrt(params.F_atm / params.F_ocean);

    void updateMomentum(const TimestepTime& tst) override
    {
#pragma omp parallel for
        for (int i = 0; i < u.rows(); ++i) {
            // Free drift ice velocity
            u(i) = uOcean(i) + NansenNumber * (uAtmos(i) * cosOceanAngle - vAtmos(i) * sinOceanAngle);
            v(i) = vOcean(i) + NansenNumber * (-uAtmos(i) * sinOceanAngle + vAtmos(i) * cosOceanAngle);
        }
    }
};
} /* namespace Nextsim */

#endif /* FREEDRIFTDYNAMICSKERNEL_HPP */
