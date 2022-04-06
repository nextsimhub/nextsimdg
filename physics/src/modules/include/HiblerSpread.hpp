/*!
 * @file HiblerSpread.hpp
 *
 * @date Apr 5, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef HIBLERSPREAD_HPP
#define HIBLERSPREAD_HPP

#include "include/Configured.hpp"
#include "ILateralIceSpread.hpp"

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

    void freeze(const TimestepTime& tsTime) override;
    void melt(const TimestepTime& tsTime) override;

private:
    static double h0;
    static double phiM;
};

}

#endif /* HIBLERSPREAD_HPP */
