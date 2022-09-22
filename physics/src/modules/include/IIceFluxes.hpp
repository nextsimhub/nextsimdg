/*!
 * @file IIceFluxes.hpp
 *
 * @date May 5, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IICEFLUXES_HPP
#define IICEFLUXES_HPP

#include "include/ModelArrayRef.hpp"
#include "include/ModelComponent.hpp"
#include "include/ModelState.hpp"
#include "include/Time.hpp"

namespace Nextsim {
//! An interface class for the calculation of ice-ocean and ice-atmosphere fluxes.
class IIceFluxes : public ModelComponent {
public:
    virtual ~IIceFluxes() = default;

    std::string getName() const override { return "IIceFluxes"; }
    void setData(const ModelState::DataMap&) override { }
    ModelState getState() const override { return ModelState(); }
    ModelState getState(const OutputLevel&) { return getState(); }

    /*!
     * Updates the ice fluxes calculation for the timestep.
     *
     * @param tStep The object containing the timestep start and duration times.
     */
    virtual void updateIce(const TimestepTime&) = 0;

protected:
    IIceFluxes()
    : qio(getSharedArray())
    , qia(getSharedArray())
    , dqia_dt(getSharedArray())
    , subl(getSharedArray())
    {
    }
    // No owned arrays
    // Shared arrays, output
    ModelArrayRef<SharedArray::Q_IO, MARBackingStore, RW> qio;
    ModelArrayRef<SharedArray::Q_IA, MARBackingStore, RW> qia;
    ModelArrayRef<SharedArray::DQIA_DT, MARBackingStore, RW> dqia_dt;
    ModelArrayRef<SharedArray::SUBLIM, MARBackingStore, RW> subl;
    // Shared arrays, input
};
}
#endif /* IICEFLUXES_HPP */
