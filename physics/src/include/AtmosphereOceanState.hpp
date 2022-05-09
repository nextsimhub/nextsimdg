/*!
 * @file AtmosphereOceanState.hpp
 *
 * @date May 9, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef PHYSICS_SRC_INCLUDE_ATMOSPHEREOCEANSTATE_HPP_
#define PHYSICS_SRC_INCLUDE_ATMOSPHEREOCEANSTATE_HPP_

#include "include/Configured.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/ModelComponent.hpp"

namespace Nextsim {

class AtmosphereOceanState : public ModelComponent, public Configured<AtmosphereOceanState> {
    AtmosphereOceanState();

    void setData(const ModelState&) override;
    ModelState getState() const override;
    ModelState getState(const OutputLevel&) const override;
    std::string getName() const override;
    std::set<std::string> hFields() const override;

    void update(const TimestepTime&);

protected:
    HField hTrueSnow;
    HField hTrueIce;

    ModelArrayRef<ProtectedArray::H_SNOW> hSnowCell;
    ModelArrayRef<ProtectedArray::H_ICE> hIceCell;
    ModelArrayRef<ProtectedArray::C_ICE> cIce;
};

} /* namespace Nextsim */

#endif /* PHYSICS_SRC_INCLUDE_ATMOSPHEREOCEANSTATE_HPP_ */
