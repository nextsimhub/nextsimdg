/*!
 * @file IOWFluxes.hpp
 *
 * @date May 5, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IOWFLUXES_HPP
#define IOWFLUXES_HPP

#include "include/ModelArrayRef.hpp"
#include "include/ModelState.hpp"
#include "include/ModelComponent.hpp"
#include "include/Time.hpp"

namespace Nextsim {
class IOWFluxes : public ModelComponent {
public:
    virtual ~IOWFluxes() = default;

    std::string getName() const override { return "IOWFluxes"; }
    void setData(const ModelState&) override { }
    ModelState getState() const override { return ModelState(); }
    ModelState getState(const OutputLevel&) { return getState(); }

    virtual void updateOW(const TimestepTime&) = 0;

protected:
    IOWFluxes() = default;
    // No owned arrays
    // Shared arrays, output
    ModelArrayRef<SharedArray::Q_OW, RW> qow;
};
}

#endif /* IOWFLUXES_HPP */
