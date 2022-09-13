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

    std::string getName() const override { return "IOceanBoundary"; }
    std::unordered_set<std::string> hFields() const override;

protected:
    HField qio; // Ice-ocean heat flux, W m⁻²
    HField sst; // Coupled or slab ocean sea surface temperature, ˚C
    HField sss; // Coupled or slab ocean sea surface salinity, PSU
    UField u; // x(east)-ward ocean current, m s⁻¹
    VField v; // y(north)-ward ocean current, m s⁻¹
};
} /* namespace Nextsim */

#endif /* IOCEANBOUNDARY_HPP */
