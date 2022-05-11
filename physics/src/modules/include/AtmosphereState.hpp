/*!
 * @file AtmosphereState.hpp
 *
 * @date May 9, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef ATMOSPHERESTATE_HPP
#define ATMOSPHERESTATE_HPP

#include "include/ModelComponent.hpp"
#include "include/Time.hpp"

namespace Nextsim {

class AtmosphereState : public ModelComponent {
public:
    AtmosphereState();
    virtual ~AtmosphereState() = default;

    void setData(const ModelState&) override;
    ModelState getState() const override;
    ModelState getState(const OutputLevel&) const override;
    std::string getName() const override;
    std::set<std::string> hFields() const override;

    void update(const TimestepTime&);

protected:
    HField tair;
    HField tdew;
    HField pair;
    HField rmix;
    HField sw_in;
    HField lw_in;
    HField snowfall;

    HField windSpeed;

    virtual void updateSpecial(const TimestepTime&) = 0;
};

} /* namespace Nextsim */

#endif /* ATMOSPHERESTATE_HPP */
