/*!
 * @file BasicIceOceanHeatFlux.hpp
 *
 * @date Oct 19, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_BASICICEOCEANHEATFLUX_HPP_
#define SRC_INCLUDE_BASICICEOCEANHEATFLUX_HPP_

#include "IIceOceanHeatFluxED.hpp"

namespace Nextsim {

//! The implementation class for the basic ice-ocean heat flux.
class BasicIceOceanHeatFluxED : public IIceOceanHeatFluxED {
public:
    BasicIceOceanHeatFluxED() = default;
    virtual ~BasicIceOceanHeatFluxED() = default;

    /*!
     * @brief Calculate the basic ice-ocean heat flux.
     *
     * @param prog PrognosticElementData for this element (constant).
     * @param exter ExternalData for this element (constant).
     * @param phys PhysicsData for this element (constant).
     * @param nsphys Nextsim physics implementation data for this element
     * (constant).
     */
    double flux(const PrognosticElementData&, const ExternalData&, const PhysicsData&,
        const NextsimPhysics&) override;
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_BASICICEOCEANHEATFLUX_HPP_ */
