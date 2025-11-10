/*!
 *
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#ifndef BASICICEOCEANHEATFLUX_HPP
#define BASICICEOCEANHEATFLUX_HPP

#include "include/IIceOceanHeatFlux.hpp"

namespace Nextsim {

//! The implementation class for the basic ice-ocean heat flux.
class BasicIceOceanHeatFlux : public IIceOceanHeatFlux {
public:
    BasicIceOceanHeatFlux()
        : IIceOceanHeatFlux()
        , mlBulkCpAccessor(getStore())
    {
    }
    virtual ~BasicIceOceanHeatFlux() = default;

    std::string getName() const override { return "BasicIceOceanHeatFlux"; }

    void update(const TimestepTime&) override;
    void updateElement(size_t i, const TimestepTime&);

protected:
    ModelArrayAccessor<Protected::ML_BULK_CP> mlBulkCpAccessor;
};

} /* namespace Nextsim */

#endif /* BASICICEOCEANHEATFLUX_HPP */
