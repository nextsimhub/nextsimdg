/*!
 * @file SlabOcean.hpp
 *
 * @date 27 Jan 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SLABOCEAN_HPP
#define SLABOCEAN_HPP

#include "include/ModelComponent.hpp"

namespace Nextsim {

class SlabOcean : public ModelComponent {
public:
    void setData(const ModelState::DataMap& ms) override;
    ModelState getState() const override;
    ModelState getState(const OutputLevel&) const override;
    std::string getName() const override { return "SlabOcean"; }

    std::unordered_set<std::string> hFields() const override;
    void update(const TimestepTime&);

private:
    HField sstSlab;
    HField sssSlab;

    static const std::string sstSlabName;
    static const std::string sssSlabName;
};

} /* namespace Nextsim */

#endif /* SLABOCEAN_HPP */
