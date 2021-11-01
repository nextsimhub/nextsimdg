/*
 * @file BasicIceOceanHeatFlux.hpp
 *
 * @date Oct 19, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_BASICICEOCEANHEATFLUX_HPP_
#define SRC_INCLUDE_BASICICEOCEANHEATFLUX_HPP_

#include "IIceOceanHeatFlux.hpp"

namespace Nextsim {

class BasicIceOceanHeatFlux : public IIceOceanHeatFlux {
public:
    BasicIceOceanHeatFlux() = default;
    virtual ~BasicIceOceanHeatFlux() = default;

    double flux(const PrognosticData&, const ExternalData&, PhysicsData&, NextsimPhysics&) override;
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_BASICICEOCEANHEATFLUX_HPP_ */
