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
#include "include/gridNames.hpp"

namespace Nextsim {
class IDynamics : public ModelComponent {
public:
    IDynamics()
        : IDynamics(false)
    {
    }
    IDynamics(bool usesDamageIn)
        : uice(ModelArray::Type::H)
        , vice(ModelArray::Type::H)
        , damage(ModelArray::Type::H)
        , hice(getStore())
        , cice(getStore())
        , hsnow(getStore())
        , damage0(getStore())
        , uwind(getStore())
        , vwind(getStore())
        , uocean(getStore())
        , vocean(getStore())
        , m_usesDamage(usesDamageIn)
    {
        getStore().registerArray(Shared::DAMAGE, &damage, RW);
    }
    virtual ~IDynamics() = default;

    ModelState getState() const override
    {
        return { {
                     { uName, mask(uice) },
                     { vName, mask(vice) },
                 },
            {} };
    }
    ModelState getState(const OutputLevel&) const override { return getState(); }
    ModelState getStateRecursive(const OutputSpec& os) const override { return os ? getState() : ModelState(); }

    std::string getName() const override { return "IDynamics"; }
    virtual void setData(const ModelState::DataMap& ms) override
    {
        uice.resize();
        vice.resize();
        damage.resize();
        if (!m_usesDamage) {
            damage = 0.;
        }
    }

    virtual void update(const TimestepTime& tst) = 0;

    /*!
     * Returns whether the dynamics implementation uses the damage field.
     */
    virtual bool usesDamage() const { return m_usesDamage; }

protected:
    // Shared ice velocity arrays
    HField uice;
    HField vice;
    // Updated damage array
    HField damage;
    // References to the DG0 finite volume data arrays
    ModelArrayRef<Shared::H_ICE, RW> hice;
    ModelArrayRef<Shared::C_ICE, RW> cice;
    ModelArrayRef<Shared::H_SNOW, RW> hsnow;
    ModelArrayRef<Protected::DAMAGE, RO> damage0;
    //ModelArrayRef<ModelComponent::SharedArray::D, MARBackingStore, RW> damage;

    // References to the forcing velocity arrays
    ModelArrayRef<Protected::WIND_U> uwind;
    ModelArrayRef<Protected::WIND_V> vwind;
    ModelArrayRef<Protected::OCEAN_U> uocean;
    ModelArrayRef<Protected::OCEAN_V> vocean;

    // Does this implementation of the dynamics use damage?
    bool m_usesDamage;
};
}

#endif /* IDYNAMICS_HPP */
