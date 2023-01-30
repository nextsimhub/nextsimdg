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
        , mlBulkCp(getProtectedArray())
    {
    }
    virtual ~BasicIceOceanHeatFlux() = default;

    void update(const TimestepTime&) override;
    void updateElement(size_t i, const TimestepTime&);
    /*!
     * Sets the relaxation time of this class.
     *
     * @param timeTIn the relaxation time in seconds to be used.
     */
    void setTimeT(double timeTIn) { timeT = timeTIn; }

protected:
    ModelArrayRef<ProtectedArray::ML_BULK_CP, MARConstBackingStore> mlBulkCp;

private:
    static double timeT; // The relaxation timescale of the sea surface temperature, s.
};

} /* namespace Nextsim */

#endif /* BASICICEOCEANHEATFLUX_HPP */
