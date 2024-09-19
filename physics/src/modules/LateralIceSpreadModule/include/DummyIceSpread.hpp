/*!
 * @file DummyIceSpread.hpp
 *
 * @date 19 Sep 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef DUMMYICESPREAD_HPP
#define DUMMYICESPREAD_HPP

#include "include/ILateralIceSpread.hpp"

namespace Nextsim {

class DummyIceSpread : public ILateralIceSpread {
public:
    DummyIceSpread()
        : ILateralIceSpread()
    {
    }
    ~DummyIceSpread() = default;

    ModelState getStateRecursive(const OutputSpec& os) const override { return ModelState(); }

    void freeze(const TimestepTime& tstep, double hice, double hsnow, double deltaHi, double newIce,
        double cice, double& qow, double& deltaCfreeze) override
    {
        deltaCfreeze = 0.;
    }
    void melt(const TimestepTime& tstep, double hice, double hsnow, double deltaHi, double cice,
        double& qow, double& deltaCmelt) override
    {
        deltaCmelt = 0.;
    }
};

} /* namespace Nextsim */

#endif /* DUMMYICESPREAD_HPP */
