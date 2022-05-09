/*
 * @file IIceOceanHeatFlux.hpp
 *
 * @date Oct 19, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_IICEOCEANHEATFLUX_HPP_
#define SRC_INCLUDE_IICEOCEANHEATFLUX_HPP_

namespace Nextsim {
class PrognosticElementData;
class ExternalData;
class PhysicsData;
class NextsimPhysics;

//! The interface class for the ice-ocean heat flux calculation.
class IIceOceanHeatFluxED {
public:
    virtual ~IIceOceanHeatFluxED() = default;

    /*!
     * @brief Calculate the ice-ocean heat flux.
     *
     * @param prog PrognosticElementData for this element (constant).
     * @param exter ExternalData for this element (constant).
     * @param phys PhysicsData for this element (constant).
     * @param nsphys Nextsim physics implementation data for this element
     * (constant).
     */
    virtual double flux(const PrognosticElementData&, const ExternalData&, const PhysicsData&,
        const NextsimPhysics&)
        = 0;
};
}
#endif /* SRC_INCLUDE_IICEOCEANHEATFLUX_HPP_ */
