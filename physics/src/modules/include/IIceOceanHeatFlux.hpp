/*
 *
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IICEOCEANHEATFLUX_HPP
#define IICEOCEANHEATFLUX_HPP

#include "include/ModelArray.hpp"
#include "include/ModelArrayAccessor.hpp"
#include "include/ModelComponent.hpp"
#include "include/Time.hpp"

namespace Nextsim {

//! The interface class for the ice-ocean heat flux calculation.
class IIceOceanHeatFlux : public ModelComponent {
public:
    IIceOceanHeatFlux()
        : sstAccessor(getStore())
        , tfAccessor(getStore())
        , ciceAccessor(getStore())
        , qioAccessor(getStore())
    {
    }
    virtual ~IIceOceanHeatFlux() = default;

    // This superclass has no state
    void setData(const ModelState::DataMap&) override {};
    std::string getName() const override { return "IIceOceanHeatFlux"; }

    /*!
     * Updates the ice ocean heat flux calculation for the timestep.
     *
     * @param tStep The object containing the timestep start and duration times.
     */
    virtual void update(const TimestepTime&) = 0;

protected:
    ModelArrayAccessor<Protected::SST> sstAccessor;
    ModelArrayAccessor<Protected::TF> tfAccessor;
    ModelArrayAccessor<Shared::C_ICE_DG> ciceAccessor;

    ModelArrayAccessor<Shared::Q_IO, RW> qioAccessor;
};
}
#endif /* IICEOCEANHEATFLUX_HPP_ */
