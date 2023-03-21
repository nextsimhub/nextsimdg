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
    IDynamics()
        : uice(ModelArray::Type::H)
        , vice(ModelArray::Type::H)
        , hice(getSharedArray())
        , cice(getSharedArray())
        , hsnow(getSharedArray())
    {
    }
    virtual ~IDynamics() = default;

    ModelState getState() const override { return ModelState(); }
    ModelState getState(const OutputLevel&) const override { return getState(); }

    std::string getName() const override { return "IDynamics"; }
    void setData(const ModelState::DataMap& ms) override
    {
        uice.resize();
        vice.resize();
    }

    virtual void update(const TimestepTime& tst) = 0;
private:
    // Shared ice velocity arrays
    HField uice;
    HField vice;
    // References to the DG0 finite volume data arrays
    ModelArrayRef<ModelComponent::SharedArray::H_ICE, MARBackingStore, RW> hice;
    ModelArrayRef<ModelComponent::SharedArray::C_ICE, MARBackingStore, RW> cice;
    ModelArrayRef<ModelComponent::SharedArray::H_SNOW, MARBackingStore, RW> hsnow;

};
}

#endif /* IDYNAMICS_HPP */
