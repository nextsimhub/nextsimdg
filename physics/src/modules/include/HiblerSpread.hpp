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

//! A class implementing the lateral spread of ice according to Hibler (1979)
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

    ModelState getStateRecursive(const OutputSpec& os) const override;

    static HelpMap& getHelpText(HelpMap& map, bool getAll);
    static HelpMap& getHelpRecursive(HelpMap&, bool getAll);

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
