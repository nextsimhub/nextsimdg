/*!
 * @file IAtmosphereBoundary.hpp
 *
 * @date Sep 12, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IATMOSPHEREBOUNDARY_HPP
#define IATMOSPHEREBOUNDARY_HPP

#include "include/ModelComponent.hpp"

namespace Nextsim {

class IAtmosphereBoundary : public ModelComponent {
    virtual ~IAtmosphereBoundary() = default;

    void setData(const ModelState::DataMap&) override;
    ModelState getState() const override;
    ModelState getState(const OutputLevel&) const override;
    ModelState getStateRecursive(const OutputSpec& os) const override;

    std::string getName() const override;
    std::unordered_set<std::string> hFields() const override;

};

} /* namespace Nextsim */

#endif /* IATMOSPHEREBOUNDARY_HPP */
