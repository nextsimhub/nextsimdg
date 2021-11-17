/*!
 * @file HiblerConcentration.hpp
 *
 * @date Nov 11, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_HIBLERCONCENTRATION_HPP
#define SRC_INCLUDE_HIBLERCONCENTRATION_HPP

#include "Configured.hpp"
#include "IConcentrationModel.hpp"

namespace Nextsim {

class HiblerConcentration : public IConcentrationModel, public Configured<HiblerConcentration> {
public:
    HiblerConcentration() = default;
    virtual ~HiblerConcentration() = default;

    void configure() override;
    enum {
        H0_KEY,
        PHIM_KEY,
    };

    double freeze(const PrognosticData&, PhysicsData&, NextsimPhysics&) const override;
    double melt(const PrognosticData&, PhysicsData&, NextsimPhysics&) const override;

    inline static void setH0(double h0_in) { h0 = h0_in; };

private:
    static double h0;
    static double phiM;
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_HIBLERCONCENTRATION_HPP */
