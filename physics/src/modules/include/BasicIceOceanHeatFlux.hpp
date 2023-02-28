/*!
 * @file BasicIceOceanHeatFlux.hpp
 *
 * @date Oct 19, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef BASICICEOCEANHEATFLUX_HPP
#define BASICICEOCEANHEATFLUX_HPP

#include "IIceOceanHeatFlux.hpp"

namespace Nextsim {

//! The implementation class for the basic ice-ocean heat flux.
class BasicIceOceanHeatFlux : public IIceOceanHeatFlux {
public:
    BasicIceOceanHeatFlux()
        : IIceOceanHeatFlux()
        , mlBulkCp(getStore())
    {
    }
    virtual ~BasicIceOceanHeatFlux() = default;

    void update(const TimestepTime&) override;
    void updateElement(size_t i, const TimestepTime&);

protected:
    ModelArrayRef<Protected::ML_BULK_CP> mlBulkCp;
};

} /* namespace Nextsim */

#endif /* BASICICEOCEANHEATFLUX_HPP */
