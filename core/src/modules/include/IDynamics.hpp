/*!
 * @file IDynamics.hpp
 *
 * @date 7 Sep 2023
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
        , hice(getStore())
        , cice(getStore())
        , hsnow(getStore())
        , uwind(getStore())
        , vwind(getStore())
        , uocean(getStore())
        , vocean(getStore())
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
    ModelArrayRef<Shared::H_ICE, RW> hice;
    ModelArrayRef<Shared::C_ICE, RW> cice;
    ModelArrayRef<Shared::H_SNOW, RW> hsnow;
    //ModelArrayRef<ModelComponent::SharedArray::D, MARBackingStore, RW> damage;

    // References to the forcing velocity arrays
    ModelArrayRef<Protected::WIND_U> uwind;
    ModelArrayRef<Protected::WIND_V> vwind;
    ModelArrayRef<Protected::OCEAN_U> uocean;
    ModelArrayRef<Protected::OCEAN_V> vocean;
};
}

#endif /* IDYNAMICS_HPP */
