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
        , uwind(getProtectedArray())
        , vwind(getProtectedArray())
        , uocean(getProtectedArray())
        , vocean(getProtectedArray())
    //, damage(getSharedArray())
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

protected:
    // Shared ice velocity arrays
    HField uice;
    HField vice;
    // References to the DG0 finite volume data arrays
    ModelArrayRef<ModelComponent::SharedArray::H_ICE, MARBackingStore, RW> hice;
    ModelArrayRef<ModelComponent::SharedArray::C_ICE, MARBackingStore, RW> cice;
    ModelArrayRef<ModelComponent::SharedArray::H_SNOW, MARBackingStore, RW> hsnow;
    // ModelArrayRef<ModelComponent::SharedArray::D, MARBackingStore, RW> damage;

    // References to the forcing velocity arrays
    ModelArrayRef<ModelComponent::ProtectedArray::WIND_U, MARConstBackingStore> uwind;
    ModelArrayRef<ModelComponent::ProtectedArray::WIND_V, MARConstBackingStore> vwind;
    ModelArrayRef<ModelComponent::ProtectedArray::OCEAN_U, MARConstBackingStore> uocean;
    ModelArrayRef<ModelComponent::ProtectedArray::OCEAN_V, MARConstBackingStore> vocean;
};
}

#endif /* IDYNAMICS_HPP */
