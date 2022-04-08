/*!
 * @file HiblerSpread.hpp
 *
 * @date Apr 5, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef HIBLERSPREAD_HPP
#define HIBLERSPREAD_HPP

#include "ILateralIceSpread.hpp"
#include "include/Configured.hpp"

namespace Nextsim {

class HiblerSpread : public ILateralIceSpread, public Configured<HiblerSpread> {
public:
    HiblerSpread()
        : ILateralIceSpread()
    {
    }
    virtual ~HiblerSpread() = default;

    void configure() override;
    enum {
        H0_KEY,
        PHIM_KEY,
    };

    //    void freeze(const TimestepTime& tsTime) override;
    //    void melt(const TimestepTime& tsTime) override;
    void freeze(const TimestepTime& tstep, double hice, double hsnow, double deltaHi, double newIce,
        double& cice, double& qow, double& deltaCfreeze) override;
    void melt(const TimestepTime& tstep, double hice, double hsnow, double deltaHi, double& cice,
        double& qow, double& deltaCmelt) override;

private:
    static double h0;
    static double phiM;
};

}

#endif /* HIBLERSPREAD_HPP */
