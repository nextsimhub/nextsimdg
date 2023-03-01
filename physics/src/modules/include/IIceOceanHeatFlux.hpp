/*
 * @file IIceOceanHeatFlux.hpp
 *
 * @date Oct 19, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IICEOCEANHEATFLUX_HPP
#define IICEOCEANHEATFLUX_HPP

#include "../../../../core/src/include/ModelArrayRef.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelComponent.hpp"
#include "include/Time.hpp"

namespace Nextsim {

//! The interface class for the ice-ocean heat flux calculation.
class IIceOceanHeatFlux : public ModelComponent {
public:
    IIceOceanHeatFlux()
        : sst(getStore())
        , tf(getStore())
        , qio(getStore())
    {
    }
    virtual ~IIceOceanHeatFlux() = default;

    // This superclass has no state
    void setData(const ModelState::DataMap&) override {};
    ModelState getState() const override { return ModelState(); }
    ModelState getState(const OutputLevel&) const override { return getState(); }
    ModelState getStateRecursive(const OutputSpec& os) const override
    {
        return os ? getState() : ModelState();
    }
    // …but it does have a name
    std::string getName() const override { return "IIceOceanHeatFlux"; }

    /*!
     * Updates the ice ocean heat flux calculation for the timestep.
     *
     * @param tStep The object containing the timestep start and duration times.
     */
    virtual void update(const TimestepTime&) = 0;

protected:
    ModelArrayRef<Protected::SST> sst;
    ModelArrayRef<Protected::TF> tf;

    ModelArrayRef<Shared::Q_IO, RW> qio;
};
}
#endif /* IICEOCEANHEATFLUX_HPP_ */
