/*!
 * @file FreeDriftDynamicsKernel.hpp
 *
 * Implementation of "classic free drift", where we ignore all \rho h terms in the momentum
 * equation. This is equivalent to assuming that the ice is very thin.
 *
 * @date 27 Mar 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Einar Ólason <einar.olason@nersc.no>
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#ifndef FREEDRIFTDYNAMICSKERNEL_HPP
#define FREEDRIFTDYNAMICSKERNEL_HPP

#include "CGDynamicsKernel.hpp"
#include "DynamicsParameters.hpp"

#include <cmath>

namespace Nextsim {

template <int DGadvection> class FreeDriftDynamicsKernel : public CGDynamicsKernel<DGadvection> {
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
    using CGDynamicsKernel<DGadvection>::updateIceOceanStress;
    using CGDynamicsKernel<DGadvection>::cosOceanAngle;
    using CGDynamicsKernel<DGadvection>::sinOceanAngle;
    using CGDynamicsKernel<DGadvection>::baseParams;

public:
    FreeDriftDynamicsKernel(const DynamicsParameters& paramsIn)
        : CGDynamicsKernel<DGadvection>(paramsIn)
        , params(paramsIn)
    {
    }

    void initialise(const ModelArray& coords, bool isSpherical, const ModelArray& mask) override
    {
        DynamicsKernel<DGadvection, DGstressComp>::initialise(coords, isSpherical, mask);

        // can't be done in the constructor since the param values are configured after construction
        FOcean = baseParams.COcean * baseParams.rhoOcean;
        FAtm = params.CAtm * params.rhoAtm;
        NansenNumber = std::sqrt(FAtm / FOcean);
    }

    void update(const TimestepTime& tst) override
    {
        updateMomentum(tst);
        applyBoundaries();

        // Let DynamicsKernel handle the advection step
        advectionAndLimits(tst);

        updateIceOceanStress(u, v);
    };

protected:
    const DynamicsParameters& params;
    double FAtm = std::numeric_limits<double>::quiet_NaN();
    double NansenNumber = std::numeric_limits<double>::quiet_NaN();
    double FOcean = std::numeric_limits<double>::quiet_NaN();

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
