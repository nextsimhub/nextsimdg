/*!
 * @file HiblerConcentration.hpp
 *
 * @date Nov 11, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_HIBLERCONCENTRATION_HPP
#define SRC_INCLUDE_HIBLERCONCENTRATION_HPP

#include "IConcentrationModel.hpp"

namespace Nextsim {

class HiblerConcentration : public IConcentrationModel {
public:
    HiblerConcentration() = default;
    virtual ~HiblerConcentration() = default;

    double freeze(const PrognosticData&, NextsimPhysics&) const override;
    double melt(const PrognosticData&, NextsimPhysics&) const override;

private:
    static double h0;
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_HIBLERCONCENTRATION_HPP */
