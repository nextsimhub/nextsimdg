/*!
 * @file IDynamics.hpp
 *
 * @date 07 Oct 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Einar Ólason <einar.olason@nersc.no>
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
        , shear(ModelArray::Type::H)
        , taux(ModelArray::Type::H)
        , tauy(ModelArray::Type::H)
        , hice(getStore())
        , cice(getStore())
        , hsnow(getStore())
        , damage0(getStore())
        , uwind(getStore())
        , vwind(getStore())
        , uocean(getStore())
        , vocean(getStore())
        , ssh(getStore())
        , m_usesDamage(usesDamageIn)
    {
        getStore().registerArray(Shared::DAMAGE, &damage, RW);
        getStore().registerArray(Protected::IO_STRESS_U, &taux, RO);
        getStore().registerArray(Protected::IO_STRESS_V, &tauy, RO);
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
    ModelState getStateRecursive(const OutputSpec& os) const override
    {
        // Ensure the base class implementation of getState() is called
        return os ? IDynamics::getState() : ModelState();
    }

    std::string getName() const override { return "IDynamics"; }
    virtual void setData(const ModelState::DataMap& ms) override
    {
        uice.resize();
        vice.resize();
        damage.resize();
        if (!m_usesDamage) {
            damage = 0.;
        }

        shear.resize();
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
    // Diagnostic outputs of shear
    HField shear;
    // Ice-ocean stress (for the coupler, mostly)
    HField taux;
    HField tauy;
    // References to the DG0 finite volume data arrays
    ModelArrayRef<Shared::H_ICE, RW> hice;
    ModelArrayRef<Shared::C_ICE, RW> cice;
    ModelArrayRef<Shared::H_SNOW, RW> hsnow;
    ModelArrayRef<Protected::DAMAGE, RO> damage0;

    // References to the forcing velocity arrays
    ModelArrayRef<Protected::WIND_U> uwind;
    ModelArrayRef<Protected::WIND_V> vwind;
    ModelArrayRef<Protected::OCEAN_U> uocean;
    ModelArrayRef<Protected::OCEAN_V> vocean;
    ModelArrayRef<Protected::SSH> ssh;

    // Does this implementation of the dynamics use damage?
    bool m_usesDamage;

    /*
     * Checks and returns if the provided data map is spherical
     */
    static bool checkSpherical(const ModelState::DataMap& ms)
    {
        // Decide between Cartesian (x & y) and spherical (longitude & latitude)
        if (ms.count(longitudeName) > 0 && ms.count(latitudeName) > 0) {
            return true;
        } else if (ms.count(xName) > 0 && ms.count(yName) > 0) {
            return false;
        } else {
            // Throw a runtime_error exception which can either be handled or not
            throw std::runtime_error("Input data must contain either Cartesian (" + xName + ", "
                + yName + ") or spherical (" + longitudeName + ", " + latitudeName
                + ") coordinates.");
        }
    }
};
}

#endif /* IDYNAMICS_HPP */
