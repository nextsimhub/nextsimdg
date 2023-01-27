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

    void update(const TimestepTime&);

private:
    HField sst;
    HField sss;
};

} /* namespace Nextsim */

#endif /* SLABOCEAN_HPP */
