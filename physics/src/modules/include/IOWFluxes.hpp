/*!
 * @file IOWFluxes.hpp
 *
 * @date May 5, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IOWFLUXES_HPP
#define IOWFLUXES_HPP

#include "include/ModelArrayRef.hpp"
#include "include/ModelComponent.hpp"
#include "include/ModelState.hpp"
#include "include/Time.hpp"

namespace Nextsim {
//! An interface class to calculate the open water ocean-atmsophere fluxes.
class IOWFluxes : public ModelComponent {
public:
    virtual ~IOWFluxes() = default;

    std::string getName() const override { return "IOWFluxes"; }
    void setData(const ModelState::DataMap&) override { }
    ModelState getState() const override { return ModelState(); }
    ModelState getState(const OutputLevel&) { return getState(); }

    /*!
     * Updates the open water fluxes for the timestep.
     *
     * @param tStep The object containing the timestep start and duration times.
     */
    virtual void updateOW(const TimestepTime&) = 0;

protected:
    IOWFluxes()
        : qow(getSharedArray())
    {
    }
    // No owned arrays
    // Shared arrays, output
    ModelArrayRef<SharedArray::Q_OW, MARBackingStore, RW> qow;
};
}

#endif /* IOWFLUXES_HPP */
