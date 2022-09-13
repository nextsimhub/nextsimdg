/*!
 * @file IOceanBoundary.hpp
 *
 * @date Sep 12, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IOCEANBOUNDARY_HPP
#define IOCEANBOUNDARY_HPP

#include "include/ModelComponent.hpp"

namespace Nextsim {

class IOceanBoundary : public ModelComponent {
    virtual ~IOceanBoundary() = default;

    void setData(const ModelState::DataMap&) override;
    ModelState getState() const override;
    ModelState getState(const OutputLevel&) const override;
    ModelState getStateRecursive(const OutputSpec& os) const override;

    std::string getName() const override;
    std::unordered_set<std::string> hFields() const override;

};
} /* namespace Nextsim */

#endif /* IOCEANBOUNDARY_HPP */
