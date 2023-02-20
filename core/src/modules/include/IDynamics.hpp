/*!
 * @file IDynamics.hpp
 *
 * @date 6 Jan 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IDYNAMICS_HPP
#define IDYNAMICS_HPP

#include "include/ModelComponent.hpp"
#include "include/Time.hpp"

namespace Nextsim {
class IDynamics : public ModelComponent {
public:
    IDynamics() = default;
    virtual ~IDynamics() = default;

    ModelState getState() const override { return ModelState(); }
    ModelState getState(const OutputLevel&) const override { return getState(); }

    std::string getName() const override { return "IDynamics"; }
    void setData(const ModelState::DataMap& ms) override = 0;

    virtual void update(const TimestepTime& tst) = 0;
};
}

#endif /* IDYNAMICS_HPP */
